!##############################################################################
!# ****************************************************************************
!# <name> list </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a double-linked list implemented as an
!# array. The list can also be used as a stack.
!#
!# The following routines are available:
!#
!# 1.) list_createList
!#     -> Create an empty list
!#
!# 2.) list_releaseList
!#     -> Release an existing list
!#
!# 3.) list_resizeList
!#     -> Reallocate memory for an existing list
!#
!# 4.) list_copyList
!#     -> Copy data to/from a list
!#
!# 5.) list_swapList
!#     -> Swap two lists
!#
!# 6.) list_getFirstInList
!#     -> Get the position of the first item in list
!#
!# 7.) list_getLastInList
!#     -> Get the position of the last item in list
!#
!# 8.) list_getNextInList
!#     -> Get the position of the next item in list
!#
!# 9.) list_prependToList = t_list_prependDble /
!#                          t_list prependSngl /
!#                          t_list_prependInt
!#     -> Prepend data to list
!#
!# 10.) list_appendToList = t_list_appendDble /
!#                          t_list_appendSngl /
!#                          t_list_appendInt
!#      -> Append data to list
!#
!# 11.) list_insertIntoList = t_list_insertDble /
!#                            t_list_insertSngl /
!#                            t_list_insertInt
!#      -> Insert data into list
!#
!# 12.) list_deleteFromList = t_list_deleteDble /
!#                            t_list_deleteSngl /
!#                            t_list_deleteInt
!#      -> Delete data from list
!#
!# 13.) list_searchInList = t_list_searchDble /
!#                          t_list_searchSngl /
!#                          t_list_searchInt
!#      -> Search for data in list
!#
!# 14.) list_printList
!#      -> Print content of list
!#
!# </purpose>
!##############################################################################

module list
  use fsystem
  use storage
  implicit none
  
  private
  public :: t_list
  public :: list_createList
  public :: list_releaseList
  public :: list_resizeList
  public :: list_copyList
  public :: list_swapList
  public :: list_getFirstInList
  public :: list_getLastInList
  public :: list_getNextInList
  public :: list_getPrevInList
  public :: list_prependToList
  public :: list_appendToList
  public :: list_insertIntoList
  public :: list_deleteFromList
  public :: list_searchInList
  public :: list_printList

!<constants>

!<constantblock description="KIND values for list data">

  ! kind value for indices in list
  integer, parameter, public :: PREC_LISTIDX = I32

!</constantblock>

!<constantblock description="Global flags for list ordering">

  ! Identifier for unordered list
  integer, parameter, public :: LIST_UNORDERED  = 0
  
  ! Identifier for increasingly ordered list
  integer, parameter, public :: LIST_INCREASING = 1

  ! Identifier for decreasingly ordered list
  integer, parameter, public :: LIST_DECREASING = 2

  ! Identifier for ordered list 
  integer, parameter, public :: LIST_CSR7       = 3

!</constantblock>

!<constantblock description="Global flags for list likn-type">

  ! Identifier fir single-linked list
  integer, parameter, public :: LIST_SINGLELINKED = 1
  
  ! Identifier for double-linked list
  integer, parameter, public :: LIST_DOUBLELINKED = 2

!</constantblock>

!<constantblock description="Global flags for list operations">

  ! Identifier for "not found in list"
  integer, parameter, public :: LIST_NOT_FOUND = -1

  ! Identifier for "found in list"
  integer, parameter, public :: LIST_FOUND     =  0
  
!</constantblock>

!<constantblock description="Internal tags for list status">
  
  ! Tag for empty list
  integer(PREC_LISTIDX), parameter, public :: LNULL =  0
  
  ! Tag for head of list
  integer(PREC_LISTIDX), parameter :: LHEAD = -2

  ! Tag for tail of list
  integer(PREC_LISTIDX), parameter :: LTAIL = -1

  ! Tag for next free item
  integer(PREC_LISTIDX), parameter :: LFREE =  0

!</constantblock>
!</constants>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<types>
!<typeblock>

  type t_list
    ! Format-Tag: Double, Single, Integer
    integer :: clistFormat        = ST_NOHANDLE
    
    ! Type-Tag: Single-linked, Double-linked
    integer :: clinkType          = LIST_UNORDERED

    ! Type of list ordering
    integer :: cordering          = LIST_UNORDERED
    
    ! Position of the last item inserted into the list
    integer :: item

    ! Number of items that are currently stored in the list
    integer(PREC_LISTIDX) :: NA   = 0

    ! Total number of items that can be stored in the list
    integer(PREC_LISTIDX) :: NNA  = 0

    ! Total number of resize operations
    integer :: NRESIZE            = 0

    ! Dimension of the auxiliary Integer values to be stored
    integer :: isizeInt            = 0

    ! Dimension of the auxiliary Double values to be stored
    integer :: isizeDble           = 0

    ! Dimension of the auxiliary Single values to be stored
    integer :: isizeSngl           = 0

    ! Factor by which the list is enlarged if new storage is allocate
    real(DP) :: dfactor           = 1.5_DP

    ! Handle to the list key
    integer :: h_Key              = ST_NOHANDLE

    ! Handle to the list next-structure
    integer :: h_Knext            = ST_NOHANDLE

    ! Handle to the list previous-structure
    integer :: h_Kprev            = ST_NOHANDLE

    ! Handle to the list auxiliary Integer data
    integer :: h_IData             = ST_NOHANDLE

    ! Handle to the list auxiliary Double data
    integer :: h_DData             = ST_NOHANDLE

    ! Handle to the list auxiliary Single data
    integer :: h_SData             = ST_NOHANDLE
    
    ! List next-structure
    ! NOTE: This array is introduced to increase performance. It
    ! should not be touched by the user. However, if the handle would
    ! be dereferenced for each operation such as search, delete,
    ! performance would be very poor.
    integer(PREC_LISTIDX), dimension(:), pointer :: Knext => null()

    ! List previous structure
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    integer(PREC_LISTIDX), dimension(:), pointer :: Kprev => null()

    ! List key data (Integer)
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    integer(PREC_LISTIDX), dimension(:), pointer :: IKey => null()

    ! List key data (Double)
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    real(DP), dimension(:), pointer :: DKey => null()

    ! List key data (Single)
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    real(SP), dimension(:), pointer :: SKey => null()
    
    ! List data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    real(DP), dimension(:,:), pointer :: DData => null()

    ! List data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    real(SP), dimension(:,:), pointer :: SData => null()

    ! List data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    integer(PREC_LISTIDX), dimension(:,:), pointer :: IData => null()
  end type t_list
  
!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  interface list_createList
    module procedure t_list_create
  end interface
  
  interface list_releaseList
    module procedure t_list_release
  end interface
  
  interface list_resizeList
    module procedure t_list_resize
  end interface
  interface resize   ! internal use
    module procedure t_list_resize
  end interface
  
  interface list_copyList
    module procedure t_list_copyto
    module procedure t_list_copytoDble
    module procedure t_list_copytoSngl
    module procedure t_list_copytoInt
    module procedure t_list_copyfrom
    module procedure t_list_copyfromDble
    module procedure t_list_copyfromSngl
    module procedure t_list_copyfromInt
  end interface
  
  interface list_swapList
    module procedure t_list_swap
  end interface
  
  interface list_getFirstInList
    module procedure t_list_first
  end interface
  
  interface list_getLastInList
    module procedure t_list_last
  end interface
  
  interface list_getNextInList
    module procedure t_list_next
  end interface

  interface list_getPrevInList
    module procedure t_list_prev
  end interface
  
  interface list_prependToList
    module procedure t_list_prependDble
    module procedure t_list_prependSngl
    module procedure t_list_prependInt
  end interface
  interface prepend   ! internal use
    module procedure t_list_prependDble
    module procedure t_list_prependSngl
    module procedure t_list_prependInt
  end interface
  
  interface list_appendToList
    module procedure t_list_appendDble
    module procedure t_list_appendSngl
    module procedure t_list_appendInt
  end interface
  interface append   ! internal use
    module procedure t_list_appendDble
    module procedure t_list_appendSngl
    module procedure t_list_appendInt
  end interface
  
  interface list_insertIntoList
    module procedure t_list_insertDble
    module procedure t_list_insertSngl
    module procedure t_list_insertInt
  end interface
  interface insert   ! internal use
    module procedure t_list_insertDble
    module procedure t_list_insertSngl
    module procedure t_list_insertInt
  end interface

  interface list_deleteFromList
    module procedure t_list_deleteDble
    module procedure t_list_deleteSngl
    module procedure t_list_deleteInt
  end interface
  interface delete   ! internal use
    module procedure t_list_deleteDble
    module procedure t_list_deleteSngl
    module procedure t_list_deleteInt
  end interface
  
  interface list_searchInList
    module procedure t_list_searchDble
    module procedure t_list_searchSngl
    module procedure t_list_searchInt
  end interface
  interface search   ! internal use
    module procedure t_list_searchDble
    module procedure t_list_searchSngl
    module procedure t_list_searchInt
  end interface
  
  interface list_printList
    module procedure t_list_print
  end interface
  
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

contains
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine t_list_create(rlist,nna,clistFormat,&
      isizeInt,isizeDble,isizeSngl,cordering,dfactor,clinkType)

!<description>
    ! This subroutine creates a new list
!</description>

!<input>
    ! Total number of items that can be stored in list
    integer(PREC_LISTIDX), intent(IN) :: nna

    ! Format-tag. Type of list format (Double,Single,Integer)
    integer, intent(IN) :: clistFormat

    ! Dimension of the auxiliary Integer values to be stored
    integer, intent(IN) :: isizeInt

    ! Dimension of the auxiliary Double values to be stored
    integer, intent(IN) :: isizeDble

    ! Dimension of the auxiliary Single values to be stored
    integer, intent(IN) :: isizeSngl

    ! OPTIONAL: Format-tag. Type of list ordering 
    integer, intent(IN), optional :: cordering

    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    real(DP), intent(IN), optional :: dfactor

    ! OPTIONAL: Type of linking (single/double). If not specified
    ! then a single-linked list is generated
    integer, intent(IN), optional :: clinkType
!</input>

!<output>
    ! list
    type(t_list), intent(OUT) :: rlist
!</output>
!</subroutine>
    
    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1_DP) rlist%dfactor=dfactor
    end if

    ! Set list format
    rlist%clistFormat=clistFormat

    ! Set ordering
    if (present(cordering)) then
      rlist%cordering=cordering
    end if

    ! Initialize list
    rlist%isizeInt  = isizeInt
    rlist%isizeSngl = isizeSngl
    rlist%isizeDble = isizeDble
    rlist%item      = LHEAD
    rlist%na        = 0
    rlist%nna       = nna
    if (present(clinkType)) then
      rlist%clinkType=clinkType
    else
      rlist%clinkType=LIST_SINGLELINKED
    end if
    
    ! Allocate memory and associate pointers
    call storage_new('t_list_create','Knext',-2,nna,ST_INT,&
        rlist%h_Knext,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rlist%h_Knext,rlist%Knext)
    
    ! Double-linked list?    
    if (rlist%clinkType == LIST_DOUBLELINKED) then
      call storage_new('t_list_create','Kprev',nna,ST_INT,&
          rlist%h_Kprev,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rlist%h_Kprev,rlist%Kprev)
    end if

    ! Allocate memory for Key
    select case(rlist%clistFormat)
    case (ST_DOUBLE)
      call storage_new('t_list_create','Key',nna,ST_DOUBLE,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double(rlist%h_Key,rlist%DKey)

    case (ST_SINGLE)
      call storage_new('t_list_create','Key',nna,ST_SINGLE,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      call storage_getbase_single(rlist%h_Key,rlist%SKey)

    case (ST_INT)
      call storage_new('t_list_create','Key',nna,ST_INT,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(rlist%h_Key,rlist%IKey)

    case DEFAULT
      print *, 't_list_create: Unsupported data format!'
      call sys_halt()
    end select
    
    ! Initialize list structure
    rlist%Knext(LFREE) = 1
    rlist%Knext(LHEAD) = LNULL
    rlist%Knext(LTAIL) = LNULL

    ! Allocate memory for auxiliary data
    if (isizeDble > 0) then
      call storage_new('t_list_create','DData',(/isizeDble,nna/),&
          ST_DOUBLE,rlist%h_DData,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2D(rlist%h_DData,rlist%DData)
    end if

     if (isizeSngl > 0) then
      call storage_new('t_list_create','SData',(/isizeSngl,nna/),&
          ST_SINGLE,rlist%h_SData,ST_NEWBLOCK_NOINIT)
      call storage_getbase_single2D(rlist%h_SData,rlist%SData)
    end if

    if (isizeInt > 0) then
      call storage_new('t_list_create','IData',(/isizeInt,nna/),&
          ST_INT,rlist%h_IData,ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2D(rlist%h_IData,rlist%IData)
    end if
  end subroutine t_list_create
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine t_list_release(rlist)

!<description>
    ! This subroutine releases an existing list
!</description>

!<inputoutput>
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>
!</subroutine>

    ! Release memory
    if (rlist%h_Key /= ST_NOHANDLE)   call storage_free(rlist%h_Key)
    if (rlist%h_Knext /= ST_NOHANDLE) call storage_free(rlist%h_Knext)
    if (rlist%h_Kprev /= ST_NOHANDLE) call storage_free(rlist%h_Kprev)

    if (rlist%h_DData /= ST_NOHANDLE) call storage_free(rlist%h_DDATA)
    if (rlist%h_SData /= ST_NOHANDLE) call storage_free(rlist%h_SDATA)
    if (rlist%h_IData /= ST_NOHANDLE) call storage_free(rlist%h_IDATA)

    nullify(rlist%Knext,rlist%Kprev,rlist%DKey,rlist%SKey,rlist%IKey)
    nullify(rlist%DData,rlist%SData,rlist%IData)

    ! Reset list
    rlist%clistFormat = ST_NOHANDLE
    rlist%clinkType   = LIST_UNORDERED
    rlist%cordering   = LIST_UNORDERED
    rlist%isizeInt    = 0
    rlist%isizeDble   = 0
    rlist%isizeSngl   = 0
    rlist%item        = 0
    rlist%na          = 0
    rlist%nna         = 0
    rlist%dfactor     = 1.5_DP
  end subroutine t_list_release

  ! ***************************************************************************

!<subroutine>
  
  subroutine t_list_resize(rlist,nna)

!<description>
    ! This subroutine reallocates memory for an existing list
!</description>

!<input>
    ! New number of total items that can be stored in the list
    integer(PREC_LISTIDX), intent(IN) :: nna
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! Set new size
    rlist%nna=nna

    ! Reallocate structures
    call storage_realloc('t_list_resize',-2,nna,rlist%h_Knext,&
        ST_NEWBLOCK_NOINIT,.true.)
    call storage_getbase_int(rlist%h_Knext,rlist%Knext)

    if (rlist%clinkType == LIST_DOUBLELINKED) then
      call storage_realloc('t_list_resize',nna,rlist%h_Kprev,&
          ST_NEWBLOCK_NOINIT,.true.)
      call storage_getbase_int(rlist%h_Kprev,rlist%Kprev)
    end if

    ! Reallocate Key
    call storage_realloc('t_list_resize',nna,rlist%h_Key,&
        ST_NEWBLOCK_NOINIT,.true.)
    select case(rlist%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(rlist%h_Key,rlist%DKey)

    case (ST_SINGLE)
      call storage_getbase_single(rlist%h_Key,rlist%SKey)

    case (ST_INT)
      call storage_getbase_int(rlist%h_Key,rlist%IKey)

    case DEFAULT
      print *, 't_list_resize: Unsupported data format!'
      call sys_halt()
    end select

    ! Reallocate auxiliary data
    if (rlist%isizeDble > 0) then
      call storage_realloc('t_list_resize',nna,rlist%h_DData,&
          ST_NEWBLOCK_NOINIT,.true.)
      call storage_getbase_double2D(rlist%h_DData,rlist%DData)
    end if

    if (rlist%isizeSngl > 0) then
      call storage_realloc('t_list_resize',nna,rlist%h_SData,&
          ST_NEWBLOCK_NOINIT,.true.)
      call storage_getbase_single2D(rlist%h_SData,rlist%SData)
    end if

    if (rlist%isizeInt > 0) then
      call storage_realloc('t_list_resize',nna,rlist%h_IData,&
          ST_NEWBLOCK_NOINIT,.true.)
      call storage_getbase_int2D(rlist%h_IData,rlist%IData)
    end if
  end subroutine t_list_resize

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_copyfrom(rlist,h_Key)

!<description>
    ! This subroutine copies the content of the list to the given
    ! handle
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist
!</input>

!<inputoutput>
    ! handle to the data
    integer, intent(INOUT) :: h_Key
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DKey
    real(SP), dimension(:), pointer :: p_SKey
    integer,  dimension(:), pointer :: p_IKey
    integer(PREC_LISTIDX) :: ipos,jpos
    
    ! Transform the content of the list to h_Key
    if (h_Key /= ST_NOHANDLE) call storage_free(h_Key)
    call storage_new('t_list_copy','Key',rlist%na,&
        rlist%clistFormat,h_Key,ST_NEWBLOCK_NOINIT)
    select case(rlist%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_Key,p_DKey)
      ipos = rlist%Knext(LHEAD)
      jpos = 0
      do
        jpos = jpos+1
        p_DKey(jpos) = rlist%DKey(ipos)
        if (ipos == rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do
      
    case (ST_SINGLE)
      call storage_getbase_single(h_Key,p_SKey)
      ipos = rlist%Knext(LHEAD)
      jpos = 0
      do
        jpos = jpos+1
        p_SKey(jpos) = rlist%SKey(ipos)
        if (ipos == rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do
      
    case (ST_INT)
      call storage_getbase_int(h_Key,p_IKey)
      ipos = rlist%Knext(LHEAD)
      jpos = 0
      do
        jpos = jpos+1
        p_IKey(jpos) = rlist%IKey(ipos)
        if (ipos == rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do
      
    case DEFAULT
      print *, 't_list_copy: Unsupported data format!'
      call sys_halt()
    end select
  end subroutine t_list_copyfrom

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_copyfromDble(rlist,p_DKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! double array
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist
!</input>

!<inputoutput>
    ! double array
    real(DP), dimension(:), intent(INOUT) :: p_DKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_LISTIDX) :: ipos,jpos

    if (rlist%clistFormat /= ST_DOUBLE) then
      print *, 't_list_copyfromDble: Unsupported data format!'
      call sys_halt()
    end if

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    do
      jpos = jpos+1
      p_DKey(jpos) = rlist%DKey(ipos)
      if (ipos == rlist%Knext(LTAIL)) exit
      ipos = rlist%Knext(ipos)
    end do
  end subroutine t_list_copyfromDble

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_copyfromSngl(rlist,p_SKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! single array
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist
!</input>

!<inputoutput>
    ! double array
    real(SP), dimension(:), intent(INOUT) :: p_SKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_LISTIDX) :: ipos,jpos

    if (rlist%clistFormat /= ST_SINGLE) then
      print *, 't_list_copyfromSngl: Unsupported data format!'
      call sys_halt()
    end if

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    do
      jpos = jpos+1
      p_SKey(jpos) = rlist%SKey(ipos)
      if (ipos == rlist%Knext(LTAIL)) exit
      ipos = rlist%Knext(ipos)
    end do
  end subroutine t_list_copyfromSngl

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_copyfromInt(rlist,p_IKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! integer array
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist
!</input>

!<inputoutput>
    ! double array
    integer, dimension(:), intent(INOUT) :: p_IKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_LISTIDX) :: ipos,jpos

    if (rlist%clistFormat /= ST_INT) then
      print *, 't_list_copyfromInt: Unsupported data format!'
      call sys_halt()
    end if

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    do
      jpos = jpos+1
      p_IKey(jpos) = rlist%IKey(ipos)
      if (ipos == rlist%Knext(LTAIL)) exit
      ipos = rlist%Knext(ipos)
    end do
  end subroutine t_list_copyfromInt

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_copyto(h_KeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given handle to the list
!</description>

!<input>
    ! handle to the data
    integer, intent(IN) :: h_KeySrc
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DKey
    real(SP), dimension(:), pointer :: p_SKey
    integer,  dimension(:), pointer :: p_IKey
    integer(PREC_LISTIDX) :: ipos,kpos
    
    ! Transform the content of h_Data to the list
    select case (rlist%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_KeySrc,p_DKey)
      do ipos=1,size(p_DKey)
        call append(rlist,p_DKey(ipos),kpos)
      end do
    case (ST_SINGLE)
      call storage_getbase_single(h_KeySrc,p_SKey)
      do ipos=1,size(p_SKey)
        call append(rlist,p_DKey(ipos),kpos)
      end do
    case (ST_INT)
      call storage_getbase_int(h_KeySrc,p_IKey)
      do ipos=1,size(p_IKey)
        call append(rlist,p_IKey(ipos),kpos)
      end do
    case DEFAULT
      print *, 't_list_copy: Unsupported data format!'
      call sys_halt()
    end select
  end subroutine t_list_copyto

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_copytoDble(p_DKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given double array to the list
!</description>

!<input>
    ! handle to the data
    real(DP), dimension(:), intent(IN) :: p_DKeySrc
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(PREC_LISTIDX) :: ipos,kpos

    if (rlist%clistFormat /= ST_DOUBLE) then
      print *, 't_list_copytoDble: Unsupported data format!'
      call sys_halt()
    end if
    
    do ipos=1,size(p_DKeySrc)
      call append(rlist,p_DKeySrc(ipos),kpos)
    end do
  end subroutine t_list_copytoDble

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_copytoSngl(p_SKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given single array to the list
!</description>

!<input>
    ! handle to the data
    real(SP), dimension(:), intent(IN) :: p_SKeySrc
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(PREC_LISTIDX) :: ipos,kpos

    if (rlist%clistFormat /= ST_SINGLE) then
      print *, 't_list_copytoSngl: Unsupported data format!'
      call sys_halt()
    end if
    
    do ipos=1,size(p_SKeySrc)
      call append(rlist,p_SKeySrc(ipos),kpos)
    end do
  end subroutine t_list_copytoSngl

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_copytoInt(p_IKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given integer array to the list
!</description>

!<input>
    ! handle to the data
    integer, dimension(:), intent(IN) :: p_IKeySrc
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(PREC_LISTIDX) :: ipos,kpos

    if (rlist%clistFormat /= ST_INT) then
      print *, 't_list_copytoInt: Unsupported data format!'
      call sys_halt()
    end if
    
    do ipos=1,size(p_IKeySrc)
      call append(rlist,p_IKeySrc(ipos),kpos)
    end do
  end subroutine t_list_copytoInt

  ! ***************************************************************************

!<subroutine>
  
  subroutine t_list_swap(rlist1,rlist2)

!<description>
    ! This subroutine swaps two lists
!</description>
    
!<inputoutput>
    ! First list
    type(t_list), intent(INOUT) :: rlist1

    ! Second list
    type(t_list), intent(INOUT) :: rlist2
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_list) :: rlist
    
    ! Check if both lists are compatible
    if (rlist1%clistFormat /= rlist2%clistFormat .or.&
        rlist1%clinkType   /= rlist2%clinkType .or. &
        rlist1%isizeInt    /= rlist2%isizeInt .or.&
        rlist1%isizeDble   /= rlist2%isizeDble .or.&
        rlist1%isizeSngl   /= rlist2%isizeSngl) then
      print *, "t_list_swap: Lists are not compatible"
      call sys_halt()
    end if

    ! Swap
    rlist  = rlist1
    rlist1 = rlist2
    rlist2 = rlist

    ! Reassociated pointers of first list
    nullify(rlist1%Knext,rlist1%Kprev,rlist1%DKey,&
        rlist1%SKey,rlist1%IKey)
    call storage_getbase_int(rlist1%h_Knext,rlist1%Knext)
    if (rlist1%h_Kprev /= ST_NOHANDLE) &
        call storage_getbase_int(rlist1%h_Kprev,rlist1%Kprev)

    select case(rlist1%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(rlist1%h_Key,rlist1%DKey)

    case (ST_SINGLE)
      call storage_getbase_single(rlist1%h_Key,rlist1%SKey)

    case (ST_INT)
      call storage_getbase_int(rlist1%h_Key,rlist1%IKey)

    case DEFAULT
      print *, 't_list_swap: Unsupported data type!'
      call sys_halt()
    end select

    if (rlist1%isizeDble > 0)&
        call storage_getbase_double2D(rlist1%h_DData,rlist1%DData)
    if (rlist1%isizeSngl > 0)&
        call storage_getbase_single2D(rlist1%h_SData,rlist1%SData)
    if (rlist1%isizeInt > 0)&
        call storage_getbase_int2D(rlist1%h_IData,rlist1%IData)
    
    ! Reassociated pointers of second list
    nullify(rlist2%Knext,rlist2%Kprev,rlist2%DKey,&
        rlist2%SKey,rlist2%IKey)
    call storage_getbase_int(rlist2%h_Knext,rlist2%Knext)
    if (rlist2%h_Kprev /= ST_NOHANDLE) &
        call storage_getbase_int(rlist2%h_Kprev,rlist2%Kprev)

    select case(rlist2%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(rlist2%h_Key,rlist2%DKey)

    case (ST_SINGLE)
      call storage_getbase_single(rlist2%h_Key,rlist2%SKey)

    case (ST_INT)
      call storage_getbase_int(rlist2%h_Key,rlist2%IKey)

    case DEFAULT
      print *, 't_list_swap: Unsupported data type!'
      call sys_halt()
    end select

    if (rlist2%isizeDble > 0)&
        call storage_getbase_double2D(rlist2%h_DData,rlist2%DData)
    if (rlist2%isizeSngl > 0)&
        call storage_getbase_single2D(rlist2%h_SData,rlist2%SData)
    if (rlist2%isizeInt > 0)&
        call storage_getbase_int2D(rlist2%h_IData,rlist2%IData)
  end subroutine t_list_swap

  ! ***************************************************************************
  
!<function>
  
  pure function t_list_first(rlist) result(ipos)

!<description>
    ! This function returns the position of the first list item
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist
!</input>

!<result>
    ! position of first item
    integer(PREC_LISTIDX) :: ipos
!</result>
!</function>
    
    ipos=rlist%Knext(LHEAD)
  end function t_list_first

  ! ***************************************************************************
  
!<function>
  
  pure function t_list_last(rlist) result(ipos)

!<description>
    ! This function returns the position of the last list item
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist
!</input>

!<result>
    ! position of last item
    integer(PREC_LISTIDX) :: ipos
!</result>
!</function>
    
    ipos=rlist%Knext(LTAIL)
  end function t_list_last

  ! ***************************************************************************

!<function>

  function t_list_next(rlist,breset) result(ipos)

!<description>
    ! This function returns the position of the next list item
    ! and resets the list if required. If the last list item is
    ! reached, then LNULL is returned and the list is reset
    ! implicitely.
!</description>

!<input>
    ! Reset list?
    logical, intent(IN) :: breset
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<result>
    ! position of last item
    integer(PREC_LISTIDX) :: ipos
!</result>
!</function>
    
    ! Reset?
    if (breset) rlist%item=LHEAD
    
    ipos = rlist%Knext(rlist%item)
    rlist%item = ipos
  end function t_list_next

   ! ***************************************************************************

!<function>

  function t_list_prev(rlist,breset) result(ipos)

!<description>
    ! This function returns the position of the previous list item
    ! and resets the list if required. This operation is only
    ! available for double-linked lists. If the first list item is
    ! reached, then LNULL is returned and the list is reset
    ! implicitely.
!</description>

!<input>
    ! Reset list?
    logical, intent(IN) :: breset
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<result>
    ! position of previous item
    integer(PREC_LISTIDX) :: ipos
!</result>
!</function>

    if (rlist%clinkType /= LIST_DOUBLELINKED) then
      print *, "t_list_prev: This operation is only available for&
          & double-linked lists!"
      call sys_halt()
    end if
    
    ! Reset?
    if (breset) rlist%item=LTAIL

    if (rlist%item == LNULL) then
      ipos = rlist%Knext(LTAIL)
    else
      ipos = rlist%Kprev(rlist%item)
    end if
    rlist%item = ipos
  end function t_list_prev

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_prependDble(rlist,dkey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends a Double data to the list
!</description>

!<input>
    ! Double key
    real(DP), intent(IN) :: dkey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat /= ST_DOUBLE) then
      print *, 't_list_prependDble: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case (rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      else
        rlist%Kprev(rlist%Knext(LHEAD)) = ipos
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
        rlist%Kprev(ipos)  = LNULL
      end if

    case DEFAULT
      print *, "t_list_prependDble: Invalid link type!"
      call sys_halt()
    end select

    ! Double key
    rlist%Dkey(ipos) = dkey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_prependDble

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_prependSngl(rlist,skey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends a Single data to the list
!</description>

!<input>
    ! Single key
    real(SP), intent(IN) :: skey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat /= ST_SINGLE) then
      print *, 't_list_prependSngl: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      else
        rlist%Kprev(rlist%Knext(LHEAD)) = ipos
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
        rlist%Kprev(ipos)  = LNULL
      end if
    case DEFAULT
      print *, "t_list_prependSngl: Invalid linktype!"
      call sys_halt()
    end select

    ! Single key
    rlist%Skey(ipos) = skey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_prependSngl

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_prependInt(rlist,ikey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends an Integer data to the list
!</description>

!<input>
    ! Integer key
    integer, intent(IN) :: ikey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat /= ST_INT) then
      print *, 't_list_prependInt: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case (rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      else
        rlist%Kprev(rlist%Knext(LHEAD)) = ipos
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
        rlist%Kprev(ipos)  = LNULL
      end if
      
    case DEFAULT
      print *, "t_list_prependInt: Invalid link type!"
      call sys_halt()
    end select

    ! Integer key
    rlist%Ikey(ipos) = ikey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_prependInt

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_appendDble(rlist,dkey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends a Double data to the list
!</description>

!<input>
    ! Double key
    real(DP), intent(IN) :: dkey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat /= ST_DOUBLE) then
      print *, 't_list_appendDble: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if
      
    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      else
        rlist%Kprev(ipos)               = rlist%Knext(LTAIL)
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if

    case DEFAULT
      print *, "t_list_appendDble: Invalid link type!"
      call sys_halt()
    end select

    ! Double key
    rlist%Dkey(ipos)   = dkey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_appendDble

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_appendSngl(rlist,skey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends a Single data to the list
!</description>

!<input>
    ! Single key
    real(SP), intent(IN) :: skey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat /= ST_SINGLE) then
      print *, 't_list_appendSngl: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      else
        rlist%Kprev(ipos)               = rlist%Knext(LTAIL)
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if

    case DEFAULT
      print *, "t_list_appendSngl: Invalid link type!"
      call sys_halt()
    end select

    ! Single key
    rlist%Skey(ipos)   = skey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_appendSngl
  
  ! ***************************************************************************

!<subroutine>

  subroutine t_list_appendInt(rlist,ikey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends an Integer data to the list
!</description>

!<input>
    ! Integer key
    integer, intent(IN) :: ikey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat /= ST_INT) then
      print *, 't_list_appendInt: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case(LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      else
        rlist%Kprev(ipos)               = rlist%Knext(LTAIL)
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if
      
    case DEFAULT
      print *, "t_list_appendInt: Invalid link type!"
      call sys_halt()
    end select

    ! Integere key
    rlist%Ikey(ipos)   = ikey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_appendInt

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_insertDble(rlist,dkey,ipred,ipos,DData,SData,IData)

!<description>
    ! This subroutine inserts a new Double data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Double key
    real(DP), intent(IN) :: dkey

    ! Position of predecessor
    integer(PREC_LISTIDX), intent(IN) :: ipred

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_DOUBLE) then
      print *, 't_list_insertDble: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      elseif (ipred == rlist%Knext(LTAIL)) then
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      elseif (ipred == rlist%Knext(LTAIL)) then
        rlist%Kprev(ipos)  = rlist%Knext(LTAIL)
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Kprev(rlist%Knext(ipred)) = ipos
        rlist%Knext(ipos)               = rlist%Knext(ipred)
        rlist%Knext(ipred)              = ipos
        rlist%Kprev(ipos)               = ipred
      end if

    case DEFAULT
      print *, "t_list_insertDble: Invalid link type!"
      call sys_halt()
    end select

    ! Double key
    rlist%Dkey(ipos)   = dkey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_insertDble

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_insertSngl(rlist,skey,ipred,ipos,DData,SData,IData)

!<description>
    ! This subroutine inserts a new Single data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Single key
    real(SP), intent(IN) :: skey

    ! Position of predecessor
    integer(PREC_LISTIDX), intent(IN) :: ipred

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_SINGLE) then
      print *, 't_list_insertSngl: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      elseif (ipred == rlist%Knext(LTAIL)) then
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      elseif (ipred == rlist%Knext(LTAIL)) then
        rlist%Kprev(ipos)  = rlist%Knext(LTAIL)
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Kprev(rlist%Knext(ipred)) = ipos
        rlist%Knext(ipos)               = rlist%Knext(ipred)
        rlist%Knext(ipred)              = ipos
        rlist%Kprev(ipos)               = ipred
      end if

    case DEFAULT
      print *, "t_list_insertSngl: Invalid link type!"
      call sys_halt()
    end select

    ! Single key
    rlist%Skey(ipos)   = skey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_insertSngl

  ! ***************************************************************************

!<subroutine>

  subroutine t_list_insertInt(rlist,ikey,ipred,ipos,DData,IData,SData)

!<description>
    ! This subroutine inserts a new Integer data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Integer key
    integer, intent(IN) :: ikey

    ! Position of predecessor
    integer(PREC_LISTIDX), intent(IN) :: ipred

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(IN), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(IN), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(IN), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer(PREC_LISTIDX), intent(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_INT) then
      print *, 't_list_insertInt: Unsupported data format!'
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call resize(rlist,ceiling(rlist%dfactor*rlist%nna))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rlist%Knext(LFREE) = ipos+1
    else
      ipos               = abs(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      elseif (ipred == rlist%Knext(LTAIL)) then
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) == LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      elseif (ipred == rlist%Knext(LTAIL)) then
        rlist%Kprev(ipos)  = rlist%Knext(LTAIL)
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Kprev(rlist%Knext(ipred)) = ipos
        rlist%Knext(ipos)               = rlist%Knext(ipred)
        rlist%Knext(ipred)              = ipos
        rlist%Kprev(ipos)               = ipred
      end if

    case DEFAULT
      print *, "t_list_insertInt: Invalid link type!"
      call sys_halt()
    end select

    ! Integer key
    rlist%Ikey(ipos)   = ikey
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) rlist%IData(:,ipos) = IData
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) rlist%DData(:,ipos) = DData
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) rlist%SData(:,ipos) = SData
  end subroutine t_list_insertInt
  
  ! ***************************************************************************
  
!<function>

  function t_list_deleteDble(rlist,dkey) result(f)

!<description>
    ! This function deletes a Double data from the list
!</description>

!<input>
    ! Data
    real(DP), intent(IN) :: dkey
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer(PREC_LISTIDX) :: ipred,ipos

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_DOUBLE) then
      print *, 't_list_deleteDble: Unsupported data format!'
      call sys_halt()
    end if

    ! Search for data
    f=search(rlist,dkey,ipred)
    if (f == LIST_NOT_FOUND) return
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    if (rlist%Knext(ipred) == rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  end function t_list_deleteDble
  
  ! ***************************************************************************
  
!<function>

  function t_list_deleteSngl(rlist,skey) result(f)

!<description>
    ! This function deletes a Single data from the list
!</description>

!<input>
    ! Data
    real(SP), intent(IN) :: skey
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer(PREC_LISTIDX) :: ipred,ipos

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_SINGLE) then
      print *, 't_list_deleteSngl: Unsupported data format!'
      call sys_halt()
    end if

    ! Search for data
    f=search(rlist,skey,ipred)
    if (f == LIST_NOT_FOUND) return
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    if (rlist%Knext(ipred) == rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  end function t_list_deleteSngl

  ! ***************************************************************************
  
!<function>

  function t_list_deleteInt(rlist,ikey) result(f)

!<description>
    ! This function deletes an Integer data from the list
!</description>

!<input>
    ! Data
    integer, intent(IN) :: ikey
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(INOUT) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer(PREC_LISTIDX) :: ipred,ipos

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_INT) then
      print *, 't_list_deleteInt: Unsupported data format!'
      call sys_halt()
    end if

    ! Search for data
    f=search(rlist,ikey,ipred)
    if (f == LIST_NOT_FOUND) return
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    if (rlist%Knext(ipred) == rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  end function t_list_deleteInt

  ! ***************************************************************************
  
!<function>

  function t_list_searchDble(rlist,dkey,ipred) result(f)

!<description>
    ! This function searches for a given Double key in the list
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist

    ! Data
    real(DP), intent(IN) :: dkey
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer(PREC_LISTIDX), intent(OUT) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer(PREC_LISTIDX) :: inext

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_DOUBLE) then
      print *, 't_list_searchDble: Unsupported data format!'
      call sys_halt()
    end if

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)

    ! Check if list is empty
    if (inext == LNULL) return

    ! What kind of ordering are we
    select case(rlist%cordering)
    case (LIST_UNORDERED)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%DKey(inext) == dkey) then
          f=LIST_FOUND; exit
        end if
        
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_INCREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%DKey(inext) == dkey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%DKey(inext) > dkey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_DECREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%DKey(inext) == dkey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%DKey(inext) < dkey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_CSR7)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%DKey(inext) == dkey) then
          f=LIST_FOUND; exit
        end if
        
        if ((ipred.ne.LHEAD) .and. rlist%DKey(inext) > dkey) exit
        
        if (rlist%Knext(ipred) == rlist%Knext(LTAIL)) then
          ipred=rlist%Knext(ipred); exit
        end if
        ipred=rlist%Knext(ipred)
      end do
    end select
  end function t_list_searchDble
  
  ! ***************************************************************************
  
!<function>

  function t_list_searchSngl(rlist,skey,ipred) result(f)

!<description>
    ! This function searches for a given Single key in the list
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist

    ! Data
    real(SP), intent(IN) :: skey
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer(PREC_LISTIDX), intent(OUT) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer(PREC_LISTIDX) :: inext

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_SINGLE) then
      print *, 't_list_searchSngl: Unsupported data format!'
      call sys_halt()
    end if

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)
    
    ! Check if list is empty
    if (inext == LNULL) return

    ! What kind of ordering are we
    select case(rlist%cordering)
    case (LIST_UNORDERED)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%SKey(inext) == skey) then
          f=LIST_FOUND; exit
        end if
        
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_INCREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%SKey(inext) == skey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%SKey(inext) > skey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_DECREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)  
        if (rlist%SKey(inext) == skey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%SKey(inext) < skey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_CSR7)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%SKey(inext) == skey) then
          f=LIST_FOUND; exit
        end if
        
        if ((ipred.ne.LHEAD) .and. rlist%SKey(inext) > skey) exit
        
        if (rlist%Knext(ipred) == rlist%Knext(LTAIL)) then
          ipred=rlist%Knext(ipred); exit
        end if
        ipred=rlist%Knext(ipred)
      end do
    end select
  end function t_list_searchSngl

  ! ***************************************************************************
  
!<function>

  function t_list_searchInt(rlist,ikey,ipred) result(f)

!<description>
    ! This function searches for a given Integer key in the list
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist

    ! Data
    integer, intent(IN) :: ikey
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer(PREC_LISTIDX), intent(OUT) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer(PREC_LISTIDX) :: inext

    ! Check if list format is ok
    if (rlist%clistFormat /= ST_INT) then
      print *, 't_list_searchInt: Unsupported data format!'
      call sys_halt()
    end if

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)
    
    ! Check if list is empty
    if (inext == LNULL) return

    ! What kind of ordering are we
    select case(rlist%cordering)
    case (LIST_UNORDERED)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%IKey(inext) == ikey) then
          f=LIST_FOUND; exit
        end if
        
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_INCREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%IKey(inext) == ikey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%IKey(inext) > ikey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_DECREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%IKey(inext) == ikey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%IKey(inext) < ikey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_CSR7)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%IKey(inext) == ikey) then
          f=LIST_FOUND; exit
        end if
        
        if ((ipred.ne.LHEAD) .and. rlist%IKey(inext) > ikey) exit
        
        if (rlist%Knext(ipred) == rlist%Knext(LTAIL)) then
          ipred=rlist%Knext(ipred); exit
        end if
        ipred=rlist%Knext(ipred)
      end do
    end select
  end function t_list_searchInt

  ! ***************************************************************************
  
!<subroutine>

  subroutine t_list_print(rlist)

!<description>
    ! This subroutine prints the content of the list
!</description>

!<input>
    ! list
    type(t_list), intent(IN) :: rlist
!</input>
!</subroutine>

    ! local variable
    integer(PREC_LISTIDX) :: ipos
    
    ipos = rlist%Knext(LHEAD)
    if (ipos == LNULL) return
    
    select case (rlist%clistFormat)
    case (ST_DOUBLE)
      do
        write(*,*) rlist%DKey(ipos)
        if (ipos == rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do

    case (ST_SINGLE)
      do
        write(*,*) rlist%SKey(ipos)
        if (ipos == rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do
      
    case (ST_INT)
      do
        write(*,*) rlist%IKey(ipos)
        if (ipos == rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do
      
    case DEFAULT
      print *, 't_list_print: Unsupported data type!'
      call sys_halt()
    end select
  end subroutine t_list_print
end module list
