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
!#  1.) list_createList
!#      -> Create an empty list
!#
!#  2.) list_releaseList
!#      -> Release an existing list
!#
!#  3.) list_resizeList
!#      -> Reallocate memory for an existing list
!#
!#  4.) list_copyToList = list_copyToListHandle /
!#                        list_copyToListDble /
!#                        list_copyToListSngl /
!#                        list_copyToListInt /
!#      -> Copy data to a list
!#
!#  5.) list_copyFromList = list_copyFromListHandle /
!#                          list_copyFromListDble /
!#                          list_copyFromListSngl /
!#                          list_copyFromListInt /
!#      -> Copy data from a list
!#
!#  6.) list_swapList
!#      -> Swap two lists
!#
!#  7.) list_getFirstInList
!#      -> Get the position of the first item in list
!#
!#  8.) list_getLastInList
!#      -> Get the position of the last item in list
!#
!#  9.) list_getNextInList
!#      -> Get the position of the next item in list
!#
!# 10.) list_prependToList = list_prependToListDble /
!#                           list_prependToListSngl /
!#                           list_prependToListInt
!#      -> Prepend data to list
!#
!# 11.) list_appendToList = list_appendToListDble /
!#                          list_appendToListSngl /
!#                          list_appendToListInt
!#      -> Append data to list
!#
!# 12.) list_insertIntoList = list_insertIntoListDble /
!#                            list_insertIntoListSngl /
!#                            list_insertIntoListInt
!#      -> Insert data into list
!#
!# 13.) list_deleteFromList = list_deleteFromListDble /
!#                            list_deleteFromListSngl /
!#                            list_deleteFromListInt
!#      -> Delete data from list
!#
!# 14.) list_searchInList = list_searchInListDble /
!#                          list_searchInListSngl /
!#                          list_searchInListInt
!#      -> Search for data in list
!#
!# 15.) list_printList
!#      -> Print content of list
!#
!# 16.) list_clearList
!#      -> Remove all content from list
!#
!# 17.) list_getByPosition = list_getByPositionDble /
!#                           list_getByPositionSngl /
!#                           list_getByPositionInt
!#
!# </purpose>
!##############################################################################

module list

  use fsystem
  use genoutput
  use storage

  implicit none
  
  private
  public :: t_list
  public :: list_createList
  public :: list_releaseList
  public :: list_resizeList
  public :: list_copyToList
  public :: list_copyFromList
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
  public :: list_clearList
  public :: list_getByPosition

!<constants>

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
  integer, parameter, public :: LNULL =  0
  
  ! Tag for head of list
  integer, parameter :: LHEAD = -2

  ! Tag for tail of list
  integer, parameter :: LTAIL = -1

  ! Tag for next free item
  integer, parameter :: LFREE =  0

!</constantblock>
!</constants>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<types>
!<typeblock>

  type t_list
    private

    ! Format-Tag: Double, Single, Integer
    integer :: clistFormat = ST_NOHANDLE
    
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

    ! Dimension of the auxiliary Integer values to be stored
    integer :: isizeInt = 0

    ! Dimension of the auxiliary Double values to be stored
    integer :: isizeDble = 0

    ! Dimension of the auxiliary Single values to be stored
    integer :: isizeSngl = 0

    ! Factor by which the list is enlarged if new storage is allocate
    real(DP) :: dfactor = 1.5_DP

    ! Handle to the list key
    integer :: h_Key = ST_NOHANDLE

    ! Handle to the list next-structure
    integer :: h_Knext = ST_NOHANDLE

    ! Handle to the list previous-structure
    integer :: h_Kprev = ST_NOHANDLE

    ! Handle to the list auxiliary Integer data
    integer :: h_IData = ST_NOHANDLE

    ! Handle to the list auxiliary Double data
    integer :: h_DData = ST_NOHANDLE

    ! Handle to the list auxiliary Single data
    integer :: h_SData = ST_NOHANDLE
    
    ! List next-structure
    ! NOTE: This array is introduced to increase performance. It
    ! should not be touched by the user. However, if the handle would
    ! be dereferenced for each operation such as search, delete,
    ! performance would be very poor.
    integer, dimension(:), pointer :: Knext => null()

    ! List previous structure
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    integer, dimension(:), pointer :: Kprev => null()

    ! List key data (Integer)
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    integer, dimension(:), pointer :: IKey => null()

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
    integer, dimension(:,:), pointer :: IData => null()
  end type t_list
  
!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************
   
  interface list_copyToList
    module procedure list_copyToListHandle
    module procedure list_copyToListDble
    module procedure list_copyToListSngl
    module procedure list_copyToListInt
  end interface

  interface list_copyFromList
    module procedure list_copyFromListHandle
    module procedure list_copyFromListDble
    module procedure list_copyFromListSngl
    module procedure list_copyFromListInt
  end interface
  
  interface list_prependToList
    module procedure list_prependToListDble
    module procedure list_prependToListSngl
    module procedure list_prependToListInt
  end interface
  
  interface list_appendToList
    module procedure list_appendToListDble
    module procedure list_appendToListSngl
    module procedure list_appendToListInt
  end interface
   
  interface list_insertIntoList
    module procedure list_insertIntoListDble
    module procedure list_insertIntoListSngl
    module procedure list_insertIntoListInt
  end interface
  
  interface list_deleteFromList
    module procedure list_deleteFromListDble
    module procedure list_deleteFromListSngl
    module procedure list_deleteFromListInt
  end interface
   
  interface list_searchInList
    module procedure list_searchInListDble
    module procedure list_searchInListSngl
    module procedure list_searchInListInt
  end interface

  interface list_getByPosition
    module procedure list_getByPositionDble
    module procedure list_getByPositionSngl
    module procedure list_getByPositionInt
  end interface
    
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

contains
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine list_createList(rlist,nna,clistFormat,&
      isizeInt,isizeDble,isizeSngl,cordering,dfactor,clinkType)

!<description>
    ! This subroutine creates a new list
!</description>

!<input>
    ! Total number of items that can be stored in list
    integer, intent(in) :: nna

    ! Format-tag. Type of list format (Double,Single,Integer)
    integer, intent(in) :: clistFormat

    ! Dimension of the auxiliary Integer values to be stored
    integer, intent(in) :: isizeInt

    ! Dimension of the auxiliary Double values to be stored
    integer, intent(in) :: isizeDble

    ! Dimension of the auxiliary Single values to be stored
    integer, intent(in) :: isizeSngl

    ! OPTIONAL: Format-tag. Type of list ordering
    integer, intent(in), optional :: cordering

    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    real(DP), intent(in), optional :: dfactor

    ! OPTIONAL: Type of linking (single/double). If not specified
    ! then a single-linked list is generated
    integer, intent(in), optional :: clinkType
!</input>

!<output>
    ! list
    type(t_list), intent(out) :: rlist
!</output>
!</subroutine>
    
    ! local variables
    integer, dimension(2) :: Isize

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
    call storage_new('list_createList','Knext',LHEAD,nna,ST_INT,&
        rlist%h_Knext,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rlist%h_Knext,rlist%Knext)
    
    ! Double-linked list?
    if (rlist%clinkType .eq. LIST_DOUBLELINKED) then
      call storage_new('list_createList','Kprev',nna,ST_INT,&
          rlist%h_Kprev,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rlist%h_Kprev,rlist%Kprev)
    end if

    ! Allocate memory for Key
    select case(rlist%clistFormat)
    case (ST_DOUBLE)
      call storage_new('list_createList','Key',nna,ST_DOUBLE,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double(rlist%h_Key,rlist%DKey)

    case (ST_SINGLE)
      call storage_new('list_createList','Key',nna,ST_SINGLE,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      call storage_getbase_single(rlist%h_Key,rlist%SKey)

    case (ST_INT)
      call storage_new('list_createList','Key',nna,ST_INT,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(rlist%h_Key,rlist%IKey)

    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_createList')
      call sys_halt()
    end select
    
    ! Initialize list structure
    rlist%Knext(LFREE) = 1
    rlist%Knext(LHEAD) = LNULL
    rlist%Knext(LTAIL) = LNULL

    ! Allocate memory for auxiliary data
    if (isizeDble > 0) then
      Isize = (/isizeDble,nna/)
      call storage_new('list_createList','DData',Isize,&
          ST_DOUBLE,rlist%h_DData,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2D(rlist%h_DData,rlist%DData)
    end if

     if (isizeSngl > 0) then
       Isize = (/isizeSngl,nna/)
      call storage_new('list_createList','SData',Isize,&
          ST_SINGLE,rlist%h_SData,ST_NEWBLOCK_NOINIT)
      call storage_getbase_single2D(rlist%h_SData,rlist%SData)
    end if

    if (isizeInt > 0) then
      Isize = (/isizeInt,nna/)
      call storage_new('list_createList','IData',Isize,&
          ST_INT,rlist%h_IData,ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2D(rlist%h_IData,rlist%IData)
    end if
  end subroutine list_createList
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine list_releaseList(rlist)

!<description>
    ! This subroutine releases an existing list
!</description>

!<inputoutput>
    type(t_list), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! Release memory
    if (rlist%h_Key .ne. ST_NOHANDLE)   call storage_free(rlist%h_Key)
    if (rlist%h_Knext .ne. ST_NOHANDLE) call storage_free(rlist%h_Knext)
    if (rlist%h_Kprev .ne. ST_NOHANDLE) call storage_free(rlist%h_Kprev)

    if (rlist%h_DData .ne. ST_NOHANDLE) call storage_free(rlist%h_DDATA)
    if (rlist%h_SData .ne. ST_NOHANDLE) call storage_free(rlist%h_SDATA)
    if (rlist%h_IData .ne. ST_NOHANDLE) call storage_free(rlist%h_IDATA)

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
  end subroutine list_releaseList

  ! ***************************************************************************

!<subroutine>
  
  subroutine list_resizeList(rlist,nna)

!<description>
    ! This subroutine reallocates memory for an existing list
!</description>

!<input>
    ! New number of total items that can be stored in the list
    integer, intent(in) :: nna
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>
!</subroutine>
    
    ! Set new size
    rlist%nna=nna

    ! Reallocate structures
    call storage_realloc('list_resizeList',LHEAD,nna,rlist%h_Knext,&
        ST_NEWBLOCK_NOINIT,.true.)
    call storage_getbase_int(rlist%h_Knext,rlist%Knext)

    if (rlist%clinkType .eq. LIST_DOUBLELINKED) then
      call storage_realloc('list_resizeList',nna,rlist%h_Kprev,&
          ST_NEWBLOCK_NOINIT,.true.)
      call storage_getbase_int(rlist%h_Kprev,rlist%Kprev)
    end if

    ! Reallocate Key
    call storage_realloc('list_resizeList',nna,rlist%h_Key,&
        ST_NEWBLOCK_NOINIT,.true.)
    select case(rlist%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(rlist%h_Key,rlist%DKey)

    case (ST_SINGLE)
      call storage_getbase_single(rlist%h_Key,rlist%SKey)

    case (ST_INT)
      call storage_getbase_int(rlist%h_Key,rlist%IKey)

    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_resizeList')
      call sys_halt()
    end select

    ! Reallocate auxiliary data
    if (rlist%isizeDble > 0) then
      call storage_realloc('list_resizeList',nna,rlist%h_DData,&
          ST_NEWBLOCK_NOINIT,.true.)
      call storage_getbase_double2D(rlist%h_DData,rlist%DData)
    end if

    if (rlist%isizeSngl > 0) then
      call storage_realloc('list_resizeList',nna,rlist%h_SData,&
          ST_NEWBLOCK_NOINIT,.true.)
      call storage_getbase_single2D(rlist%h_SData,rlist%SData)
    end if

    if (rlist%isizeInt > 0) then
      call storage_realloc('list_resizeList',nna,rlist%h_IData,&
          ST_NEWBLOCK_NOINIT,.true.)
      call storage_getbase_int2D(rlist%h_IData,rlist%IData)
    end if
  end subroutine list_resizeList

  ! ***************************************************************************

!<subroutine>

  subroutine list_copyFromListHandle(rlist,h_Key)

!<description>
    ! This subroutine copies the content of the list to the given
    ! handle
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist
!</input>

!<inputoutput>
    ! handle to the data
    integer, intent(inout) :: h_Key
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DKey
    real(SP), dimension(:), pointer :: p_SKey
    integer,  dimension(:), pointer :: p_IKey
    integer :: ipos,jpos
    
    ! Transform the content of the list to h_Key
    if (h_Key .ne. ST_NOHANDLE) call storage_free(h_Key)
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
        if (ipos .eq. rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do
      
    case (ST_SINGLE)
      call storage_getbase_single(h_Key,p_SKey)
      ipos = rlist%Knext(LHEAD)
      jpos = 0
      do
        jpos = jpos+1
        p_SKey(jpos) = rlist%SKey(ipos)
        if (ipos .eq. rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do
      
    case (ST_INT)
      call storage_getbase_int(h_Key,p_IKey)
      ipos = rlist%Knext(LHEAD)
      jpos = 0
      do
        jpos = jpos+1
        p_IKey(jpos) = rlist%IKey(ipos)
        if (ipos .eq. rlist%Knext(LTAIL)) exit
        ipos = rlist%Knext(ipos)
      end do
      
    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_copyFromListHandle')
      call sys_halt()
    end select
  end subroutine list_copyFromListHandle

  ! ***************************************************************************

!<subroutine>

  subroutine list_copyFromListDble(rlist,p_DKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! double array
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist
!</input>

!<inputoutput>
    ! double array
    real(DP), dimension(:), intent(inout) :: p_DKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos,jpos

    if (rlist%clistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_copyFromListDble')
      call sys_halt()
    end if

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    do
      jpos = jpos+1
      p_DKey(jpos) = rlist%DKey(ipos)
      if (ipos .eq. rlist%Knext(LTAIL)) exit
      ipos = rlist%Knext(ipos)
    end do
  end subroutine list_copyFromListDble

  ! ***************************************************************************

!<subroutine>

  subroutine list_copyFromListSngl(rlist,p_SKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! single array
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist
!</input>

!<inputoutput>
    ! double array
    real(SP), dimension(:), intent(inout) :: p_SKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos,jpos

    if (rlist%clistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_copyFromListSngl')
      call sys_halt()
    end if

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    do
      jpos = jpos+1
      p_SKey(jpos) = rlist%SKey(ipos)
      if (ipos .eq. rlist%Knext(LTAIL)) exit
      ipos = rlist%Knext(ipos)
    end do
  end subroutine list_copyFromListSngl

  ! ***************************************************************************

!<subroutine>

  subroutine list_copyFromListInt(rlist,p_IKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! integer array
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist
!</input>

!<inputoutput>
    ! double array
    integer, dimension(:), intent(inout) :: p_IKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos,jpos

    if (rlist%clistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_copyFromListInt')
      call sys_halt()
    end if

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    do
      jpos = jpos+1
      p_IKey(jpos) = rlist%IKey(ipos)
      if (ipos .eq. rlist%Knext(LTAIL)) exit
      ipos = rlist%Knext(ipos)
    end do
  end subroutine list_copyFromListInt

  ! ***************************************************************************

!<subroutine>

  subroutine list_copyToListHandle(h_KeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given handle to the list
!</description>

!<input>
    ! handle to the data
    integer, intent(in) :: h_KeySrc
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DKey
    real(SP), dimension(:), pointer :: p_SKey
    integer,  dimension(:), pointer :: p_IKey
    integer :: ipos,kpos
    
    ! Transform the content of h_Data to the list
    select case (rlist%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_KeySrc,p_DKey)
      do ipos=1,size(p_DKey)
        call list_appendToList(rlist,p_DKey(ipos),kpos)
      end do
    case (ST_SINGLE)
      call storage_getbase_single(h_KeySrc,p_SKey)
      do ipos=1,size(p_SKey)
        call list_appendToList(rlist,p_DKey(ipos),kpos)
      end do
    case (ST_INT)
      call storage_getbase_int(h_KeySrc,p_IKey)
      do ipos=1,size(p_IKey)
        call list_appendToList(rlist,p_IKey(ipos),kpos)
      end do
    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_copyToListHandle')
      call sys_halt()
    end select
  end subroutine list_copyToListHandle

  ! ***************************************************************************

!<subroutine>

  subroutine list_copyToListDble(p_DKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given double array to the list
!</description>

!<input>
    ! handle to the data
    real(DP), dimension(:), intent(in) :: p_DKeySrc
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ipos,kpos

    if (rlist%clistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_copyToListDble')
      call sys_halt()
    end if
    
    do ipos=1,size(p_DKeySrc)
      call list_appendToList(rlist,p_DKeySrc(ipos),kpos)
    end do
  end subroutine list_copyToListDble

  ! ***************************************************************************

!<subroutine>

  subroutine list_copyToListSngl(p_SKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given single array to the list
!</description>

!<input>
    ! handle to the data
    real(SP), dimension(:), intent(in) :: p_SKeySrc
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ipos,kpos

    if (rlist%clistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_copyToListSngl')
      call sys_halt()
    end if
    
    do ipos=1,size(p_SKeySrc)
      call list_appendToList(rlist,p_SKeySrc(ipos),kpos)
    end do
  end subroutine list_copyToListSngl

  ! ***************************************************************************

!<subroutine>

  subroutine list_copyToListInt(p_IKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given integer array to the list
!</description>

!<input>
    ! handle to the data
    integer, dimension(:), intent(in) :: p_IKeySrc
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ipos,kpos

    if (rlist%clistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_copyToListInt')
      call sys_halt()
    end if
    
    do ipos=1,size(p_IKeySrc)
      call list_appendToList(rlist,p_IKeySrc(ipos),kpos)
    end do
  end subroutine list_copyToListInt

  ! ***************************************************************************

!<subroutine>
  
  subroutine list_swapList(rlist1,rlist2)

!<description>
    ! This subroutine swaps two lists
!</description>
    
!<inputoutput>
    ! First list
    type(t_list), intent(inout) :: rlist1

    ! Second list
    type(t_list), intent(inout) :: rlist2
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_list) :: rlist
    
    ! Check if both lists are compatible
    if (rlist1%clistFormat .ne. rlist2%clistFormat .or.&
        rlist1%clinkType   .ne. rlist2%clinkType .or. &
        rlist1%isizeInt    .ne. rlist2%isizeInt .or.&
        rlist1%isizeDble   .ne. rlist2%isizeDble .or.&
        rlist1%isizeSngl   .ne. rlist2%isizeSngl) then

      call output_line('Lists are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_swapList')
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
    if (rlist1%h_Kprev .ne. ST_NOHANDLE) &
        call storage_getbase_int(rlist1%h_Kprev,rlist1%Kprev)

    select case(rlist1%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(rlist1%h_Key,rlist1%DKey)

    case (ST_SINGLE)
      call storage_getbase_single(rlist1%h_Key,rlist1%SKey)

    case (ST_INT)
      call storage_getbase_int(rlist1%h_Key,rlist1%IKey)

    case DEFAULT
      call output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_swapList')
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
    if (rlist2%h_Kprev .ne. ST_NOHANDLE) &
        call storage_getbase_int(rlist2%h_Kprev,rlist2%Kprev)

    select case(rlist2%clistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(rlist2%h_Key,rlist2%DKey)

    case (ST_SINGLE)
      call storage_getbase_single(rlist2%h_Key,rlist2%SKey)

    case (ST_INT)
      call storage_getbase_int(rlist2%h_Key,rlist2%IKey)

    case DEFAULT
      call output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_swapList')
      call sys_halt()
    end select

    if (rlist2%isizeDble > 0)&
        call storage_getbase_double2D(rlist2%h_DData,rlist2%DData)
    if (rlist2%isizeSngl > 0)&
        call storage_getbase_single2D(rlist2%h_SData,rlist2%SData)
    if (rlist2%isizeInt > 0)&
        call storage_getbase_int2D(rlist2%h_IData,rlist2%IData)
  end subroutine list_swapList

  ! ***************************************************************************
  
!<function>
  
  pure function list_getFirstInList(rlist) result(ipos)

!<description>
    ! This function returns the position of the first list item
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist
!</input>

!<result>
    ! position of first item
    integer :: ipos
!</result>
!</function>
    
    ipos=rlist%Knext(LHEAD)
  end function list_getFirstInList

  ! ***************************************************************************
  
!<function>
  
  pure function list_getLastInList(rlist) result(ipos)

!<description>
    ! This function returns the position of the last list item
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist
!</input>

!<result>
    ! position of last item
    integer :: ipos
!</result>
!</function>
    
    ipos=rlist%Knext(LTAIL)
  end function list_getLastInList

  ! ***************************************************************************

!<function>

  function list_getNextInList(rlist,breset) result(ipos)

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
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! position of last item
    integer :: ipos
!</result>
!</function>
    
    ! Reset?
    if (breset) rlist%item=LHEAD
    
    ipos = rlist%Knext(rlist%item)
    rlist%item = ipos
  end function list_getNextInList

   ! ***************************************************************************

!<function>

  function list_getPrevInList(rlist,breset) result(ipos)

!<description>
    ! This function returns the position of the previous list item
    ! and resets the list if required. This operation is only
    ! available for double-linked lists. If the first list item is
    ! reached, then LNULL is returned and the list is reset
    ! implicitely.
!</description>

!<input>
    ! Reset list?
    logical, intent(in) :: breset
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! position of previous item
    integer :: ipos
!</result>
!</function>

    if (rlist%clinkType .ne. LIST_DOUBLELINKED) then
      call output_line('This operation is only available for double-linked lists!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_getPrevInList')
      call sys_halt()
    end if
    
    ! Reset?
    if (breset) rlist%item=LTAIL

    if (rlist%item .eq. LNULL) then
      ipos = rlist%Knext(LTAIL)
    else
      ipos = rlist%Kprev(rlist%item)
    end if
    rlist%item = ipos
  end function list_getPrevInList

  ! ***************************************************************************

!<subroutine>

  subroutine list_prependToListDble(rlist,dkey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends a Double data to the list
!</description>

!<input>
    ! Double key
    real(DP), intent(in) :: dkey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_prependToListDble')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_prependToListDble')
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
  end subroutine list_prependToListDble

  ! ***************************************************************************

!<subroutine>

  subroutine list_prependToListSngl(rlist,skey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends a Single data to the list
!</description>

!<input>
    ! Single key
    real(SP), intent(in) :: skey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_prependToListSngl')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_prependToListSngl')
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
  end subroutine list_prependToListSngl

  ! ***************************************************************************

!<subroutine>

  subroutine list_prependToListInt(rlist,ikey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends an Integer data to the list
!</description>

!<input>
    ! Integer key
    integer, intent(in) :: ikey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_prependToListInt')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_prependToListInt')
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
  end subroutine list_prependToListInt

  ! ***************************************************************************

!<subroutine>

  subroutine list_appendToListDble(rlist,dkey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends a Double data to the list
!</description>

!<input>
    ! Double key
    real(DP), intent(in) :: dkey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_appendToListDble')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if
      
    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_appendToListDble')
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
  end subroutine list_appendToListDble

  ! ***************************************************************************

!<subroutine>

  subroutine list_appendToListSngl(rlist,skey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends a Single data to the list
!</description>

!<input>
    ! Single key
    real(SP), intent(in) :: skey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_appendToListSngl')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_appendToListSngl')
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
  end subroutine list_appendToListSngl
  
  ! ***************************************************************************

!<subroutine>

  subroutine list_appendToListInt(rlist,ikey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends an Integer data to the list
!</description>

!<input>
    ! Integer key
    integer, intent(in) :: ikey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_appendToListInt')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_appendToListInt')
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
  end subroutine list_appendToListInt

  ! ***************************************************************************

!<subroutine>

  subroutine list_insertIntoListDble(rlist,dkey,ipred,ipos,DData,SData,IData)

!<description>
    ! This subroutine inserts a new Double data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Double key
    real(DP), intent(in) :: dkey

    ! Position of predecessor
    integer, intent(in) :: ipred

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_insertIntoListDble')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      elseif (ipred .eq. rlist%Knext(LTAIL)) then
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      elseif (ipred .eq. rlist%Knext(LTAIL)) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_insertIntoListDble')
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
  end subroutine list_insertIntoListDble

  ! ***************************************************************************

!<subroutine>

  subroutine list_insertIntoListSngl(rlist,skey,ipred,ipos,DData,SData,IData)

!<description>
    ! This subroutine inserts a new Single data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Single key
    real(SP), intent(in) :: skey

    ! Position of predecessor
    integer, intent(in) :: ipred

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_insertIntoListSngl')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      elseif (ipred .eq. rlist%Knext(LTAIL)) then
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      elseif (ipred .eq. rlist%Knext(LTAIL)) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_insertIntoListSngl')
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
  end subroutine list_insertIntoListSngl

  ! ***************************************************************************

!<subroutine>

  subroutine list_insertIntoListInt(rlist,ikey,ipred,ipos,DData,IData,SData)

!<description>
    ! This subroutine inserts a new Integer data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Integer key
    integer, intent(in) :: ikey

    ! Position of predecessor
    integer, intent(in) :: ipred

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_insertIntoListInt')
      call sys_halt()
    end if
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resizeList(rlist,ceiling(rlist%dfactor*rlist%nna))
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
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      elseif (ipred .eq. rlist%Knext(LTAIL)) then
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      else
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%Knext(LHEAD) .eq. LNULL) then
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      elseif (ipred .eq. rlist%Knext(LTAIL)) then
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
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_insertIntoListInt')
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
  end subroutine list_insertIntoListInt
  
  ! ***************************************************************************
  
!<function>

  function list_deleteFromListDble(rlist,dkey) result(f)

!<description>
    ! This function deletes a Double data from the list
!</description>

!<input>
    ! Data
    real(DP), intent(in) :: dkey
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_deleteFromListDble')
      call sys_halt()
    end if

    ! Search for data
    f=list_searchInList(rlist,dkey,ipred)
    if (f .eq. LIST_NOT_FOUND) return
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    if (rlist%Knext(ipred) .eq. rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  end function list_deleteFromListDble
  
  ! ***************************************************************************
  
!<function>

  function list_deleteFromListSngl(rlist,skey) result(f)

!<description>
    ! This function deletes a Single data from the list
!</description>

!<input>
    ! Data
    real(SP), intent(in) :: skey
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_deleteFromListSngl')
      call sys_halt()
    end if

    ! Search for data
    f=list_searchInList(rlist,skey,ipred)
    if (f .eq. LIST_NOT_FOUND) return
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    if (rlist%Knext(ipred) .eq. rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  end function list_deleteFromListSngl

  ! ***************************************************************************
  
!<function>

  function list_deleteFromListInt(rlist,ikey) result(f)

!<description>
    ! This function deletes an Integer data from the list
!</description>

!<input>
    ! Data
    integer, intent(in) :: ikey
!</input>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_deleteFromListInt')
      call sys_halt()
    end if

    ! Search for data
    f=list_searchInList(rlist,ikey,ipred)
    if (f .eq. LIST_NOT_FOUND) return
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    if (rlist%Knext(ipred) .eq. rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  end function list_deleteFromListInt

  ! ***************************************************************************
  
!<function>

  function list_searchInListDble(rlist,dkey,ipred) result(f)

!<description>
    ! This function searches for a given Double key in the list
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist

    ! Data
    real(DP), intent(in) :: dkey
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer, intent(out) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: inext

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_searchInListDble')
      call sys_halt()
    end if

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)

    ! Check if list is empty
    if (inext .eq. LNULL) return

    ! What kind of ordering are we
    select case(rlist%cordering)
    case (LIST_UNORDERED)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%DKey(inext) .eq. dkey) then
          f=LIST_FOUND; exit
        end if
        
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_INCREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%DKey(inext) .eq. dkey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%DKey(inext) > dkey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_DECREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%DKey(inext) .eq. dkey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%DKey(inext) < dkey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_CSR7)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%DKey(inext) .eq. dkey) then
          f=LIST_FOUND; exit
        end if
        
        if ((ipred.ne.LHEAD) .and. rlist%DKey(inext) > dkey) exit
        
        if (rlist%Knext(ipred) .eq. rlist%Knext(LTAIL)) then
          ipred=rlist%Knext(ipred); exit
        end if
        ipred=rlist%Knext(ipred)
      end do
    end select
  end function list_searchInListDble
  
  ! ***************************************************************************
  
!<function>

  function list_searchInListSngl(rlist,skey,ipred) result(f)

!<description>
    ! This function searches for a given Single key in the list
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist

    ! Data
    real(SP), intent(in) :: skey
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer, intent(out) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: inext

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_searchInListSngl')
      call sys_halt()
    end if

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)
    
    ! Check if list is empty
    if (inext .eq. LNULL) return

    ! What kind of ordering are we
    select case(rlist%cordering)
    case (LIST_UNORDERED)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%SKey(inext) .eq. skey) then
          f=LIST_FOUND; exit
        end if
        
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_INCREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%SKey(inext) .eq. skey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%SKey(inext) > skey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_DECREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%SKey(inext) .eq. skey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%SKey(inext) < skey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_CSR7)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%SKey(inext) .eq. skey) then
          f=LIST_FOUND; exit
        end if
        
        if ((ipred.ne.LHEAD) .and. rlist%SKey(inext) > skey) exit
        
        if (rlist%Knext(ipred) .eq. rlist%Knext(LTAIL)) then
          ipred=rlist%Knext(ipred); exit
        end if
        ipred=rlist%Knext(ipred)
      end do
    end select
  end function list_searchInListSngl

  ! ***************************************************************************
  
!<function>

  function list_searchInListInt(rlist,ikey,ipred) result(f)

!<description>
    ! This function searches for a given Integer key in the list
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist

    ! Data
    integer, intent(in) :: ikey
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer, intent(out) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: inext

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_searchInListInt')
      call sys_halt()
    end if

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)
    
    ! Check if list is empty
    if (inext .eq. LNULL) return

    ! What kind of ordering are we
    select case(rlist%cordering)
    case (LIST_UNORDERED)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%IKey(inext) .eq. ikey) then
          f=LIST_FOUND; exit
        end if
        
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_INCREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%IKey(inext) .eq. ikey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%IKey(inext) > ikey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_DECREASING)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%IKey(inext) .eq. ikey) then
          f=LIST_FOUND; exit
        end if
        
        if (rlist%IKey(inext) < ikey) exit
        ipred=rlist%Knext(ipred)
      end do
      
    case (LIST_CSR7)
      do while(ipred.ne.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        if (rlist%IKey(inext) .eq. ikey) then
          f=LIST_FOUND; exit
        end if
        
        if ((ipred.ne.LHEAD) .and. rlist%IKey(inext) > ikey) exit
        
        if (rlist%Knext(ipred) .eq. rlist%Knext(LTAIL)) then
          ipred=rlist%Knext(ipred); exit
        end if
        ipred=rlist%Knext(ipred)
      end do
    end select
  end function list_searchInListInt

  ! ***************************************************************************
  
!<subroutine>

  subroutine list_printList(rlist)

!<description>
    ! This subroutine prints the content of the list
!</description>

!<input>
    ! list
    type(t_list), intent(in) :: rlist
!</input>
!</subroutine>

    ! local variable
    integer :: ipos
    
    ipos = rlist%Knext(LHEAD)
    if (ipos .eq. LNULL) return
    
    select case (rlist%clistFormat)
    case (ST_DOUBLE)
      do while (ipos .ne. rlist%Knext(LTAIL))
        write(*,*) rlist%DKey(ipos)
        ipos = rlist%Knext(ipos)
      end do

    case (ST_SINGLE)
      do while (ipos .ne. rlist%Knext(LTAIL))
        write(*,*) rlist%SKey(ipos)
        ipos = rlist%Knext(ipos)
      end do
      
    case (ST_INT)
      do while (ipos .ne. rlist%Knext(LTAIL))
        write(*,*) rlist%IKey(ipos)
        ipos = rlist%Knext(ipos)
      end do
      
    case DEFAULT
      call output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_printList')
      call sys_halt()
    end select
  end subroutine list_printList

  ! ***************************************************************************
  
!<subroutine>

  subroutine list_clearList(rlist)

!<description>
    ! This subroutine clears the content of the list
!</description>

!<inputoutput>
    ! list
    type(t_list), intent(inout) :: rlist
!</input>
!</subroutine>

    ! Initialize list structure
    rlist%item         = LHEAD
    rlist%na           = 0
    rlist%Knext(LFREE) = 1
    rlist%Knext(LHEAD) = LNULL
    rlist%Knext(LTAIL) = LNULL
    
  end subroutine list_clearList

  ! ***************************************************************************

!<subroutine>

  subroutine list_getByPositionDble(rlist, ipos, dkey, DData, SData, IData)

!<description>
    ! This subroutine returns the key and data stored at the given position
!</description>

!<input>
    ! List
    type(t_list), intent(in) :: rlist

    ! Position of the data
    integer, intent(in) :: ipos
!</input>

!<output>
    ! Double key
    real(DP), intent(out) :: dkey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(out), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(out), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(out), optional :: IData
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_getByPositionDble')
      call sys_halt()
    end if

    ! Check if position is valid
    if (ipos > rlist%na) then
      call output_line('Invalid position!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_getByPositionDble')
      call sys_halt()
    end if

    ! Double key
    dkey = rlist%Dkey(ipos)
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) IData = rlist%IData(:,ipos)
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) DData = rlist%DData(:,ipos)
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) SData = rlist%SData(:,ipos)

  end subroutine list_getByPositionDble

  ! ***************************************************************************

!<subroutine>

  subroutine list_getByPositionSngl(rlist, ipos, skey, DData, SData, IData)

!<description>
    ! This subroutine returns the key and data stored at the given position
!</description>

!<input>
    ! List
    type(t_list), intent(in) :: rlist

    ! Position of the data
    integer, intent(in) :: ipos
!</input>

!<output>
    ! Single key
    real(SP), intent(out) :: skey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(out), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(out), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(out), optional :: IData
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_getByPositionSngl')
      call sys_halt()
    end if

    ! Check if position is valid
    if (ipos > rlist%na) then
      call output_line('Invalid position!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_getByPositionSngl')
      call sys_halt()
    end if

    ! Double key
    skey = rlist%Skey(ipos)
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) IData = rlist%IData(:,ipos)
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) DData = rlist%DData(:,ipos)
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) SData = rlist%SData(:,ipos)

  end subroutine list_getByPositionSngl

  ! ***************************************************************************

!<subroutine>

  subroutine list_getByPositionInt(rlist, ipos, ikey, DData, SData, IData)

!<description>
    ! This subroutine returns the key and data stored at the given position
!</description>

!<input>
    ! List
    type(t_list), intent(in) :: rlist

    ! Position of the data
    integer, intent(in) :: ipos
!</input>

!<output>
    ! Integer key
    integer, intent(out) :: ikey

    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(out), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(out), optional :: SData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(out), optional :: IData
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rlist%clistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_getByPositionInt')
      call sys_halt()
    end if

    ! Check if position is valid
    if (ipos > rlist%na) then
      call output_line('Invalid position!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_getByPositionInt')
      call sys_halt()
    end if

    ! Double key
    ikey = rlist%Ikey(ipos)
    
    ! Optional data
    if ((rlist%isizeInt > 0) .and. &
        present(IData)) IData = rlist%IData(:,ipos)
    
    if ((rlist%isizeDble > 0) .and. &
        present(DData)) DData = rlist%DData(:,ipos)
    
    if ((rlist%isizeSngl > 0) .and. &
        present(SData)) SData = rlist%SData(:,ipos)

  end subroutine list_getByPositionInt

end module list
