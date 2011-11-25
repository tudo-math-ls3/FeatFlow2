!##############################################################################
!# ****************************************************************************
!# <name> arraylist </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements an array of (linear) linked lists, also
!# known as arraylist, which is implemented as one large storage block
!#
!# The following routines are available:
!#
!# 1.) arrlst_createArrayList = arrlst_createArrayList /
!#                              arrlst_createArrayListTbl
!#     -> Create an empty arraylist
!#
!# 2.) arrlst_releaseArrayList = arrlst_releaseArrayList /
!#                               arrlst_releaseArrayListTbl
!#     -> Release an existing arraylist
!#
!# 3.) arrlst_resizeArrayList
!#     -> Reallocate memory for an existing arraylist
!#
!# 4.) arrlst_copyArrayList = arrlst_copyFromArrayList /
!#                            arrlst_copyFromArrayListDble /
!#                            arrlst_copyFromArrayListSngl /
!#                            arrlst_copyFromArrayListInt /
!#                            arrlst_copyToArrayList /
!#                            arrlst_copyToArrayListDble /
!#                            arrlst_copyToArrayListSngl /
!#                            arrlst_copyToArrayListInt
!#     -> Copy data to/from an arraylist for a given table
!#
!# 5.) arrlst_copyArrayListTable = arrlst_copyFromArrayListTbl /
!#                                 arrlst_copyFromArrayListDbleTbl /
!#                                 arrlst_copyFromArrayListSnglTbl /
!#                                 arrlst_copyFromArrayListIntTbl /
!#                                 arrlst_copyToArrayListTbl /
!#                                 arrlst_copyToArrayListDbleTbl /
!#                                 arrlst_copyToArrayListSnglTbl /
!#                                 arrlst_copyToArrayListIntTbl
!#     -> Copy data to/from an arraylist for a complete table
!#        structure
!#
!# 6.) arrlst_swapArrayList
!#     -> Swap two lists in the table
!#
!# 7.) arrlst_getFirstInArrayList
!#     -> Get the position of the first item in arraylist
!#
!# 8.) arrlst_getLastInArrayList
!#     -> Get the position of the last item in arraylist
!#
!# 9.) arrlst_getNextInArrayList
!#     -> Get the position of the next item in arraylist
!#
!# 10.) arrlst_prependToArrayList = arrlst_prependToArrayListDble /
!#                                  arrlst_prependToArrayListSngl /
!#                                  arrlst_prependToArrayListInt
!#     -> Prepend data to arraylist
!#
!# 11.) arrlst_appendToArrayist = arrlst_appendToArrayListDble /
!#                                arrlst_appendToArrayListSngl /
!#                                arrlst_appendToArrayListInt
!#      -> Append data to arraylist
!#
!# 12.) arrlst_insertIntoArrayList = arrlst_insertIntoArrayListDble /
!#                                   arrlst_insertIntoArrayListSngl /
!#                                   arrlst_insertIntoArrayListInt
!#      -> Insert data into arraylist
!#
!# 13.) arrlst_deleteFromArrayList = arrlst_deleteFromArrayListDble /
!#                                   arrlst_deleteFromArrayListSngl /
!#                                   arrlst_deleteFromArrayListInt
!#      -> Delete data from arraylist
!#
!# 14.) arrlst_searchInArrayList = arrlst_searchInArrayListDble /
!#                                 arrlst_searchInArrayListSngl /
!#                                 arrlst_searchInArrayListInt
!#      -> Search for data in arraylist
!#
!# 15.) arrlst_printArrayList
!#      -> Print content of arraylist
!#
!# 16.) arrlst_infoArrayList
!#      -> Print information about arraylist
!#
!# 17.) arrlst_duplicateArrayList
!#      -> Create a duplicate / backup of an arraylist.
!#
!# 18.) arrlst_restoreArrayList
!#      -> Restore an arraylist from a previous backup
!# </purpose>
!##############################################################################

module arraylist
  use fsystem
  use storage
  use genoutput
  implicit none
  
  private
  public :: t_arraylist
  public :: arrlst_createArrayList
  public :: arrlst_releaseArrayList
  public :: arrlst_resizeArrayList
  public :: arrlst_copyArrayList
  public :: arrlst_copyArrayListTable
  public :: arrlst_swapArrayList
  public :: arrlst_getFirstInArrayList
  public :: arrlst_getLastInArrayList
  public :: arrlst_getNextInArrayList
  public :: arrlst_prependToArrayList
  public :: arrlst_appendToArrayList
  public :: arrlst_insertIntoArrayList
  public :: arrlst_deleteFromArrayList
  public :: arrlst_searchInArrayList
  public :: arrlst_printArrayList
  public :: arrlst_infoArrayList
  public :: arrlst_duplicateArrayList
  public :: arrlst_restoreArrayList

!<constants>

!<constantblock description="Global flags for arraylist ordering">

  ! Identifier for unordered arraylist
  integer, parameter, public :: ARRAYLIST_UNORDERED  = 0
  
  ! Identifier for increasingly ordered list
  integer, parameter, public :: ARRAYLIST_INCREASING = 1

  ! Identifier for decreasingly ordered arraylist
  integer, parameter, public :: ARRAYLIST_DECREASING = 2

  ! Identifier for ordered arraylist
  integer, parameter, public :: ARRAYLIST_CSR7       = 3

!</constantblock>

!<constantblock description="Global flags for arraylist operations">

  ! Identifier for "not found in arraylist"
  integer, parameter, public :: ARRAYLIST_NOT_FOUND = -1

  ! Identifier for "found in arraylist"
  integer, parameter, public :: ARRAYLIST_FOUND     =  0
  
!</constantblock>

!<constantblock description="Internal tags for arraylist status">
  
  ! Tag for empty arraylist
  integer, parameter, public :: ARRLST_NULL          =  0

  ! Tag for next free position in storage of arraylist
  integer, parameter :: ARRLST_FREE                  = 0

  ! Tag for head of each list
  integer, parameter :: ARRLST_HEAD                  = 1

  ! Tag for tail of each list
  integer, parameter :: ARRLST_TAIL                  = 2

  ! Tag for last item stored in each list
  integer, parameter :: ARRLST_ITEM                  = 3

  ! Tag for number of entries stored in each list
  integer, parameter :: ARRLST_NA                    = 4

!</constantblock>
!</constants>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<types>
!<typeblock>

  type t_arraylist
    ! Format-Tag:
    integer :: carraylistFormat                   = ST_NOHANDLE
    
    ! Type of arraylist ordering
    integer :: cordering                          = ARRAYLIST_UNORDERED
    
    ! Number of tables that are currently stored in the arraylist
    integer :: NTABLE              = 0
    
    ! Total number of tables that can be stored in the arraylist
    integer :: NNTABLE             = 0

    ! Total number of tables that can initially be stored in the
    ! arraylist. This information is needed to compute the growth of
    ! the arraylist after several resize operations
    integer :: NNTABLE0            = 0

    ! Number of items that are currently stored in all lists
    integer :: NA              = 0

    ! Total number of items that can be stored in all lists
    integer :: NNA             = 0

    ! Total number of items that can initially be stored in all
    ! lists. This information is needed to compute the growth of the
    ! arraylist after several resize operations
    integer :: NNA0            = 0

    ! Number of resize operations performed with the arraylist
    integer :: NRESIZE                            = 0

    ! Factor by which the arraylist is enlarged if new storage is allocate
    real(DP) :: dfactor                           = 1.5_DP

    ! Handle to the lookup table
    integer :: h_Ktable                           = ST_NOHANDLE

    ! Handle to the arraylist structure
    integer :: h_Knext                            = ST_NOHANDLE
    
    ! Handle to the arraylist data
    integer :: h_Data                             = ST_NOHANDLE

    ! Table structure
    ! NOTE: This array is introduced to increase performance. It
    ! should not be touched by the user. However, if the handle would
    ! be dereferenced for each operation such as search, delete,
    ! performance would be very poor.
    integer, dimension(:,:), pointer :: p_Ktable => null()
    
    ! ArrayList structure
    ! NOTE: This array is introduced to increase performance (see above).
    integer, dimension(:), pointer ::   p_Knext => null()

    ! ArrayList data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    real(DP), dimension(:), pointer ::                     p_DData => null()

    ! ArrayList data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    real(SP), dimension(:), pointer ::                     p_FData => null()

    ! ArrayList data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    integer, dimension(:), pointer ::   p_IData => null()
  end type t_arraylist
  
!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************
  interface arrlst_createArrayList
    module procedure arrlst_createArrayListDefault
    module procedure arrlst_createArrayListTbl
  end interface

  interface arrlst_releaseArrayList
    module procedure arrlst_releaseArrayListDefault
    module procedure arrlst_releaseArrayListTbl
  end interface
  
  interface arrlst_copyArrayList
    module procedure arrlst_copyFromArrayList
    module procedure arrlst_copyFromArrayListDble
    module procedure arrlst_copyFromArrayListSngl
    module procedure arrlst_copyFromArrayListInt
    module procedure arrlst_copyToArrayList
    module procedure arrlst_copyToArrayListDble
    module procedure arrlst_copyToArrayListSngl
    module procedure arrlst_copyToArrayListInt
  end interface

  interface arrlst_copyArrayListTable
    module procedure arrlst_copyFromArrayListTbl
    module procedure arrlst_copyFromArrayListDbleTbl
    module procedure arrlst_copyFromArrayListSnglTbl
    module procedure arrlst_copyFromArrayListIntTbl
    module procedure arrlst_copyToArrayListTbl
    module procedure arrlst_copyToArrayListDbleTbl
    module procedure arrlst_copyToArrayListSnglTbl
    module procedure arrlst_copyToArrayListIntTbl
  end interface
    
  interface arrlst_prependToArrayList
    module procedure arrlst_prependToArrayListDble
    module procedure arrlst_prependToArrayListSngl
    module procedure arrlst_prependToArrayListInt
  end interface
   
  interface arrlst_appendToArrayList
    module procedure arrlst_appendToArrayListDble
    module procedure arrlst_appendToArrayListSngl
    module procedure arrlst_appendToArrayListInt
  end interface
   
  interface arrlst_insertIntoArrayList
    module procedure arrlst_insertIntoArrayListDble
    module procedure arrlst_insertIntoArrayListSngl
    module procedure arrlst_insertIntoArrayListInt
  end interface
 
  interface arrlst_deleteFromArrayList
    module procedure arrlst_deleteFromArrayListDble
    module procedure arrlst_deleteFromArrayListSngl
    module procedure arrlst_deleteFromArrayListInt
  end interface
    
  interface arrlst_searchInArrayList
    module procedure arrlst_searchInArrayListDble
    module procedure arrlst_searchInArrayListSngl
    module procedure arrlst_searchInArrayListInt
  end interface
  
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

contains
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine arrlst_createArrayListDefault(rarraylist,nntable,nna,carraylistFormat,cordering,dfactor)

!<description>
    ! This subroutine creates a new arraylist
!</description>

!<input>
    ! Total number of tables
    integer, intent(in) :: nntable

    ! Total number of items that can be stored in all lists
    integer, intent(in) :: nna

    ! Format-tag. Type of list format (Double,Single,Integer)
    integer, intent(in) :: carraylistFormat

    ! OPTIONAL: Format-tag. Type of list ordering
    integer, intent(in), optional :: cordering
    
    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    real(DP), intent(in), optional :: dfactor
!</input>

!<output>
    ! list
    type(t_arraylist), intent(out) :: rarraylist
!</output>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    
    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1_DP) rarraylist%dfactor=dfactor
    end if

    ! Set list format
    rarraylist%carraylistFormat=carraylistFormat

    ! Set ordering
    if (present(cordering)) then
      rarraylist%cordering=cordering
    end if

    ! Initialize list
    rarraylist%NTABLE   = 0
    rarraylist%NNTABLE  = nntable
    rarraylist%NNTABLE0 = nntable
    rarraylist%NA       = 0
    rarraylist%NNA      = nna
    rarraylist%NNA0     = nna
    rarraylist%NRESIZE  = 0

    ! Allocate memory and associate pointers
    Isize=(/4,nntable/)
    call storage_new('arrlst_createArrayListDefault','Ktable',Isize,&
        ST_INT,rarraylist%h_Ktable,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int2D(rarraylist%h_Ktable,rarraylist%p_Ktable)
    
    call storage_new('arrlst_createArrayListDefault','Knext',ARRLST_FREE,nna,ST_INT,&
        rarraylist%h_Knext,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rarraylist%h_Knext,rarraylist%p_Knext)

    select case(rarraylist%carraylistFormat)
    case (ST_DOUBLE)
      call storage_new('arrlst_createArrayListDefault','Data',nna,ST_DOUBLE,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double(rarraylist%h_Data,rarraylist%p_DData)
      
    case (ST_SINGLE)
      call storage_new('arrlst_createArrayListDefault','Data',nna,ST_SINGLE,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      call storage_getbase_single(rarraylist%h_Data,rarraylist%p_FData)
      
    case (ST_INT)
      call storage_new('arrlst_createArrayListDefault','Data',nna,ST_INT,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(rarraylist%h_Data,rarraylist%p_IData)
      
    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_createArrayListDefault')
      call sys_halt()
    end select
    
    ! Initialize list structures
    rarraylist%p_Knext(ARRLST_FREE)    = 1
  end subroutine arrlst_createArrayListDefault
  
  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_createArrayListTbl(rarraylist,itable)

!<description>
    ! This subroutine creates new table entries up to position
    ! itable. Note: This subroutine cannot be called from outside of
    ! this module by the user. It is only used internally if some
    ! item should be added to a list for which there exists no table.
!</description>

!<input>
    ! Number of last table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! Check if table number is valid
    if (itable < 1) then
      call output_line('Invalid table number',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_createArrayListTbl')
      call sys_halt()
    end if

    ! Resize tables if required
    if (itable > rarraylist%NNTABLE) call arrlst_resizeArrayListTbl(&
        rarraylist,ceiling(itable*rarraylist%dfactor))
    
    ! Initialize structures
    rarraylist%p_Ktable(ARRLST_HEAD:ARRLST_NA,rarraylist%NTABLE+1:itable) = ARRLST_NULL

    ! Set new table size
    rarraylist%NTABLE = max(rarraylist%NTABLE,itable)
  end subroutine arrlst_createArrayListTbl

  ! ***************************************************************************

!<subroutine>
  
  subroutine arrlst_releaseArrayListDefault(rarraylist)

!<description>
    ! This subroutine releases an existing arraylist
!</description>

!<inputoutput>
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Release memory
    if (rarraylist%h_Ktable .ne. ST_NOHANDLE) call storage_free(rarraylist%h_Ktable)
    if (rarraylist%h_Knext .ne. ST_NOHANDLE)  call storage_free(rarraylist%h_Knext)
    if (rarraylist%h_Data .ne. ST_NOHANDLE)   call storage_free(rarraylist%h_Data)
    nullify(rarraylist%p_Ktable,rarraylist%p_Knext,rarraylist%p_DData,&
        rarraylist%p_FData,rarraylist%p_IData)

    ! Reset list
    rarraylist%carraylistFormat = ST_NOHANDLE
    rarraylist%cordering   = ARRAYLIST_UNORDERED
    rarraylist%NTABLE      = 0
    rarraylist%NNTABLE     = 0
    rarraylist%NNTABLE0    = 0
    rarraylist%NA          = 0
    rarraylist%NNA         = 0
    rarraylist%NNA0        = 0
    rarraylist%dfactor     = 1.5_DP
    rarraylist%NRESIZE     = 0
  end subroutine arrlst_releaseArrayListDefault

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_releaseArrayListTbl(rarraylist,itable)

!<description>
    ! This subroutine releases a table from the arraylist
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Check if table exists
    if (itable < 1 .or. itable > rarraylist%NTABLE) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_releaseArrayListTbl')
      call sys_halt()
    end if

    ! Decrease number of entries by the number of entries present
    ! in the table which is released
    rarraylist%NA = rarraylist%NA - rarraylist%p_Ktable(ARRLST_NA,itable)

    ! Reset table
    rarraylist%p_Ktable(ARRLST_HEAD:ARRLST_NA,itable) = ARRLST_NULL

    ! Decrease number of tables if the last table has been deleted
    if (itable .eq. rarraylist%NTABLE)&
        rarraylist%NTABLE = rarraylist%NTABLE-1
  end subroutine arrlst_releaseArrayListTbl

  ! ***************************************************************************

!<subroutine>
  
  subroutine arrlst_resizeArrayList(rarraylist,nna)

!<description>
    ! This subroutine reallocates memory for an existing list
!</description>

!<input>
    ! New number of total items that can be stored in the list
    integer, intent(in) :: nna
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! Set new size and increase counter
    rarraylist%NNA     = nna
    rarraylist%NRESIZE = rarraylist%NRESIZE+1

    call storage_realloc('arrlst_resizeArrayList',nna+1,&
        rarraylist%h_Knext,ST_NEWBLOCK_NOINIT,.true.)
    call storage_realloc('arrlst_resizeArrayList',nna,&
        rarraylist%h_Data,ST_NEWBLOCK_NOINIT,.true.)
    call storage_getbase_int(rarraylist%h_Knext,rarraylist%p_Knext)

    select case(rarraylist%carraylistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(rarraylist%h_Data,rarraylist%p_DData)

    case (ST_SINGLE)
      call storage_getbase_single(rarraylist%h_Data,rarraylist%p_FData)

    case (ST_INT)
      call storage_getbase_int(rarraylist%h_Data,rarraylist%p_IData)

    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_resizeArrayList')
      call sys_halt()
    end select
  end subroutine arrlst_resizeArrayList

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_resizeArrayListTbl(rarraylist,nntable)

!<description>
    ! This subroutine reallocates memory for the lookup table
!</description>

!<input>
    ! New number of tables
    integer, intent(in) :: nntable
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Set new size
    rarraylist%NNTABLE = nntable
    rarraylist%NRESIZE = rarraylist%NRESIZE+1

    call storage_realloc('arrlst_resizeArrayListTbl',nntable,&
        rarraylist%h_Ktable,ST_NEWBLOCK_NOINIT,.true.)
    call storage_getbase_int2D(rarraylist%h_Ktable,rarraylist%p_Ktable)

  end subroutine arrlst_resizeArrayListTbl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyFromArrayList(rarraylist,itable,h_Data)

!<description>
    ! This subroutine copies the content of the list of a given table
    ! to the given handle.
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! Number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! handle to the data
    integer, intent(inout) :: h_Data
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DData
    real(SP), dimension(:), pointer :: p_FData
    integer,  dimension(:), pointer :: p_IData
    integer :: isize

    ! Transform the content of the list to h_Data
    if (h_Data .eq. ST_NOHANDLE) then
      call storage_new('arrlst_copyFromArrayList','Data',rarraylist%NA,&
          rarraylist%carraylistFormat,h_Data,ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_Data,isize)
      if (isize < rarraylist%NA) then
        call storage_realloc('arrlst_copyFromArrayList',rarraylist%NA,&
            h_Data,ST_NEWBLOCK_NOINIT,.false.)
      end if
    end if
    
    ! What kind of data are we?
    select case(rarraylist%carraylistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_Data,p_DData)
      call arrlst_copyArrayList(rarraylist,itable,p_DData)
      
    case (ST_SINGLE)
      call storage_getbase_single(h_Data,p_FData)
      call arrlst_copyArrayList(rarraylist,itable,p_FData)
      
    case (ST_INT)
      call storage_getbase_int(h_Data,p_IData)
      call arrlst_copyArrayList(rarraylist,itable,p_IData)
      
    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayList')
      call sys_halt()
    end select
  end subroutine arrlst_copyFromArrayList

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyFromArrayListTbl(rarraylist,h_Data,h_Table)

!<description>
    ! This subroutine copies the content of the table and all lists
    ! to the given handles.
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist

!</input>

!<inputoutput>
    ! handle to the data
    integer, intent(inout) :: h_Data

    ! handle to the table
    integer, intent(inout) :: h_Table
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Table
    real(DP), dimension(:), pointer :: p_DData
    real(SP), dimension(:), pointer :: p_FData
    integer,  dimension(:), pointer :: p_IData
    integer :: isize

    ! Transform the content of the arraylist and the
    ! table to h_Data and h_Table, respectively
    if (h_Table .eq. ST_NOHANDLE) then
      call storage_new('arrlst_copyFromArrayListTbl','Table',&
          rarraylist%NTABLE+1,ST_INT,h_Table,ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_Table,isize)
      if (isize < rarraylist%NTABLE+1) then
        call storage_realloc('arrlst_copyFromArrayListTbl',&
            rarraylist%NTABLE+1,h_Table,ST_NEWBLOCK_NOINIT,.false.)
      end if
    end if
    call storage_getbase_int(h_Table,p_Table)
    
    if (h_Data .eq. ST_NOHANDLE) then
      call storage_new('arrlst_copyFromArrayListTbl','Data',&
          rarraylist%NA,rarraylist%carraylistFormat,h_Data,ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_Data,isize)
      if (isize < rarraylist%NA) then
        call storage_realloc('arrlst_copyFromArrayListTbl',&
            rarraylist%NA,h_Data,ST_NEWBLOCK_NOINIT,.false.)
      end if
    end if

    ! What kind of array are we?
    select case(rarraylist%carraylistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_Data,p_DData)
      call arrlst_copyArrayListTable(rarraylist,p_DData,p_Table)

    case (ST_SINGLE)
      call storage_getbase_single(h_Data,p_FData)
      call arrlst_copyArrayListTable(rarraylist,p_FData,p_Table)

    case (ST_INT)
      call storage_getbase_int(h_Data,p_IData)
      call arrlst_copyArrayListTable(rarraylist,p_IData,p_Table)
      
    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListTbl')
      call sys_halt()
    end select
  end subroutine arrlst_copyFromArrayListTbl
  
  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyFromArrayListDble(rarraylist,itable,p_DData,ndata)

!<description>
    ! This subroutine copies the content of the list of the given
    ! table to the given double array
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! Number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! double array
    real(DP), dimension(:), intent(inout) :: p_DData
!</inputoutput>

!<output>
    ! OPTIONAL: number of data entries
    integer, intent(out), optional :: ndata
!</output>
!</subroutine>

    ! local variables
    integer :: ipos,icount

    ! Check if table is valid
    if (itable < 1 .or. itable > rarraylist%NTABLE) then
      call output_line('Invalid table number',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListDble')
      call sys_halt()
    end if
    
    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListDble')
      call sys_halt()
    end if

    icount = 0
    ipos = rarraylist%p_Ktable(ARRLST_HEAD,itable)
    do
      icount = icount+1
      p_DData(icount) = rarraylist%p_DData(ipos)
      if (ipos .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) exit
      ipos = rarraylist%p_Knext(ipos)
    end do

    if (present(ndata)) ndata=icount
  end subroutine arrlst_copyFromArrayListDble

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyFromArrayListDbleTbl(rarraylist,p_DData,p_Table)

!<description>
    ! This subroutine copies the content of the table and all lists
    ! to the given arrays.
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<inputoutput>
    ! double array
    real(DP), dimension(:), intent(inout) :: p_DData

    ! table array
    integer, dimension(:), intent(inout) :: p_Table
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: icount,itable,ntable,ipos

    ! Check if table array is valid
    ntable=size(p_Table)-1
    if (ntable+1 < rarraylist%NTABLE) then
      call output_line('Invalid dimension of table array!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListDbleTbl')
      call sys_halt()
    end if
    
    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListDbleTbl')
      call sys_halt()
    end if

    icount=1
    do itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%p_Ktable(ARRLST_HEAD,itable)
      do while (ipos .ne. ARRLST_NULL)
        p_DData(icount) = rarraylist%p_DData(ipos)
        icount = icount+1
        if (ipos .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) exit
        ipos = rarraylist%p_Knext(ipos)
      end do
    end do
    p_Table(ntable+1)=icount+1
  end subroutine arrlst_copyFromArrayListDbleTbl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyFromArrayListSngl(rarraylist,itable,p_FData,ndata)

!<description>
    ! This subroutine copies the content of the list of the given
    ! table to the given single array
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist
    
    ! Number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! double array
    real(SP), dimension(:), intent(inout) :: p_FData
!</inputoutput>

!<output>
    ! OPTIONAL: number of data entries
    integer, intent(out), optional :: ndata
!</output>
!</subroutine>

    ! local variables
    integer :: ipos,icount

    ! Check if table is valid
    if (itable < 1 .or. itable > rarraylist%NTABLE) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListSngl')
      call sys_halt()
    end if

    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListSngl')
      call sys_halt()
    end if

    icount = 0
    ipos = rarraylist%p_Ktable(ARRLST_HEAD,itable)
    do
      icount = icount+1
      p_FData(icount) = rarraylist%p_FData(ipos)
      if (ipos .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) exit
      ipos = rarraylist%p_Knext(ipos)
    end do

    if (present(ndata)) ndata=icount
  end subroutine arrlst_copyFromArrayListSngl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyFromArrayListSnglTbl(rarraylist,p_FData,p_Table)

!<description>
    ! This subroutine copies the content of the table and all lists
    ! to the given arrays.
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<inputoutput>
    ! single array
    real(SP), dimension(:), intent(inout) :: p_FData

    ! table array
    integer, dimension(:), intent(inout) :: p_Table
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: icount,itable,ntable,ipos

    ! Check if table array is valid
    ntable=size(p_Table)-1
    if (ntable+1 < rarraylist%NTABLE) then
      call output_line('Invalid dimension of table array!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListSnglTbl')
      call sys_halt()
    end if
    
    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListSnglTbl')
      call sys_halt()
    end if

    icount=1
    do itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%p_Ktable(ARRLST_HEAD,itable)
      do while(ipos .ne. ARRLST_NULL)
        p_FData(icount) = rarraylist%p_FData(ipos)
        icount = icount+1
        if (ipos .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) exit
        ipos = rarraylist%p_Knext(ipos)
      end do
    end do
    p_Table(ntable+1)=icount+1
  end subroutine arrlst_copyFromArrayListSnglTbl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyFromArrayListInt(rarraylist,itable,p_IData,ndata)

!<description>
    ! This subroutine copies the content of the list of the given
    ! table to the given integer array.
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! Number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! double array
    integer, dimension(:), intent(inout) :: p_IData
!</inputoutput>

!<output>
    ! OPTIONAL: number of data entries
    integer, intent(out), optional :: ndata
!</output>
!</subroutine>

    ! local variables
    integer :: ipos,icount

    ! Check if table is valid
    if (itable < 1 .or. itable > rarraylist%NTABLE) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListInt')
      call sys_halt()
    end if

    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListInt')
      call sys_halt()
    end if

    icount = 0
    ipos = rarraylist%p_Ktable(ARRLST_HEAD,itable)
    do
      icount = icount+1
      p_IData(icount) = rarraylist%p_IData(ipos)
      if (ipos .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) exit
      ipos = rarraylist%p_Knext(ipos)
    end do

    if (present(ndata)) ndata=icount
  end subroutine arrlst_copyFromArrayListInt

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyFromArrayListIntTbl(rarraylist,p_IData,p_Table)

!<description>
    ! This subroutine copies the content of the table and all lists
    ! to the given arrays.
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<inputoutput>
    ! integer array
    integer, dimension(:), intent(inout) :: p_IData

    ! table array
    integer, dimension(:), intent(inout) :: p_Table
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: icount,itable,ntable,ipos

    ! Check if table array is valid
    ntable=size(p_Table)-1
    if (ntable+1 < rarraylist%NTABLE) then
      call output_line('Invalid dimension of table array!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListIntTbl')
      call sys_halt()
    end if
    
    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListIntTbl')
      call sys_halt()
    end if

    icount=1
    do itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%p_Ktable(ARRLST_HEAD,itable)
      do while(ipos .ne. ARRLST_NULL)
        p_IData(icount) = rarraylist%p_IData(ipos)
        icount=icount+1
        if (ipos .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) exit
        ipos = rarraylist%p_Knext(ipos)
      end do
    end do
    p_Table(ntable+1)=icount
  end subroutine arrlst_copyFromArrayListIntTbl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyToArrayList(h_DataSrc,itable,rarraylist)

!<description>
    ! This subroutine copies the content of the given handle to the
    ! list associated with the given table
!</description>

!<input>
    ! handle to the data
    integer, intent(in) :: h_DataSrc
     
    ! Number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DData
    real(SP), dimension(:), pointer :: p_FData
    integer,  dimension(:), pointer :: p_IData
    
    ! Transform the content of h_Data to the list
    select case (rarraylist%carraylistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_DataSrc,p_DData)
      call arrlst_copyArrayList(p_DData,itable,rarraylist)

    case (ST_SINGLE)
      call storage_getbase_single(h_DataSrc,p_FData)
      call arrlst_copyArrayList(p_FData,itable,rarraylist)

    case (ST_INT)
      call storage_getbase_int(h_DataSrc,p_IData)
      call arrlst_copyArrayList(p_IData,itable,rarraylist)
      stop

    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayList')
      call sys_halt()
    end select
  end subroutine arrlst_copyToArrayList

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyToArrayListTbl(h_DataSrc,rarraylist,h_Table)

!<description>
    ! This subroutine copies the content of the given handle to the
    ! lists of the arraylist making use of the second handle to
    ! generate the table
!</description>

!<input>
    ! handle to the data
    integer, intent(in) :: h_DataSrc

    ! handle to the table
    integer, intent(in) :: h_Table
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    integer,  dimension(:), pointer :: p_Table
    real(DP), dimension(:), pointer :: p_DData
    real(SP), dimension(:), pointer :: p_FData
    integer,  dimension(:), pointer :: p_IData
    
    ! Set pointer to table
    call storage_getbase_int(h_Table,p_Table)

    ! Transform the content of h_Data to the list
    select case (rarraylist%carraylistFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_DataSrc,p_DData)
      call arrlst_copyArrayListTable(p_DData,rarraylist,p_Table)

    case (ST_SINGLE)
      call storage_getbase_single(h_DataSrc,p_FData)
      call arrlst_copyArrayListTable(p_FData,rarraylist,p_Table)

    case (ST_INT)
      call storage_getbase_int(h_DataSrc,p_IData)
      call arrlst_copyArrayListTable(p_IData,rarraylist,p_Table)

    case DEFAULT
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListTbl')
      call sys_halt()
    end select
  end subroutine arrlst_copyToArrayListTbl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyToArrayListDble(p_DDataSrc,itable,rarraylist)

!<description>
    ! This subroutine copies the content of the given double array to
    ! the list associated with the given table
!</description>

!<input>
    ! pointer to the data
    real(DP), dimension(:), intent(in) :: p_DDataSrc

    ! number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ipos,kpos

    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListDble')
      call sys_halt()
    end if
    
    do ipos=1,size(p_DDataSrc)
      call arrlst_appendToArrayList(rarraylist,itable,p_DDataSrc(ipos),kpos)
    end do
  end subroutine arrlst_copyToArrayListDble

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyToArrayListDbleTbl(p_DDataSrc,rarraylist,p_Table)

!<description>
    ! This subroutine copies the content of the given double array to
    ! the lists of the arraylist making use of the second table array
!</description>

!<input>
    ! pointer to the data
    real(DP), dimension(:), intent(in) :: p_DDataSrc

    ! pointer to the table
    integer, dimension(:), intent(in) :: p_Table
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: itable,ntable,ipos,kpos

    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListDbleTbl')
      call sys_halt()
    end if
    
    ntable=size(p_Table)-1
    do itable=1,ntable
      do ipos=p_Table(itable),p_Table(itable+1)-1
        call arrlst_appendToArrayList(rarraylist,itable,p_DDataSrc(ipos),kpos)
      end do
    end do
  end subroutine arrlst_copyToArrayListDbleTbl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyToArrayListSngl(p_FDataSrc,itable,rarraylist)

!<description>
    ! This subroutine copies the content of the given single array to
    ! the list associated with the given table
!</description>

!<input>
    ! pointer to the data
    real(SP), dimension(:), intent(in) :: p_FDataSrc

    ! number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ipos,kpos

    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListSngl')
      call sys_halt()
    end if
    
    do ipos=1,size(p_FDataSrc)
      call arrlst_appendToArrayList(rarraylist,itable,p_FDataSrc(ipos),kpos)
    end do
  end subroutine arrlst_copyToArrayListSngl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyToArrayListSnglTbl(p_FDataSrc,rarraylist,p_Table)

!<description>
    ! This subroutine copies the content of the given single array to
    ! the lists of the arraylist making use of the second table array
!</description>

!<input>
    ! pointer to the data
    real(SP), dimension(:), intent(in) :: p_FDataSrc

    ! pointer to the table
    integer, dimension(:), intent(in) :: p_Table
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: itable,ntable,ipos,kpos

    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListSnglTbl')
      call sys_halt()
    end if

    ntable=size(p_Table)-1
    do itable=1,ntable
      do ipos=p_Table(itable),p_Table(itable+1)-1
        call arrlst_appendToArrayList(rarraylist,itable,p_FDataSrc(ipos),kpos)
      end do
    end do
  end subroutine arrlst_copyToArrayListSnglTbl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyToArrayListInt(p_IDataSrc,itable,rarraylist)

!<description>
    ! This subroutine copies the content of the given integer array
    ! to the list associated with the given table
!</description>

!<input>
    ! pointer to the data
    integer, dimension(:), intent(in) :: p_IDataSrc

    ! number of table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ipos,kpos

    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListInt')
      call sys_halt()
    end if
    
    do ipos=1,size(p_IDataSrc)
      call arrlst_appendToArrayList(rarraylist,itable,p_IDataSrc(ipos),kpos)
    end do
  end subroutine arrlst_copyToArrayListInt

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_copyToArrayListIntTbl(p_IDataSrc,rarraylist,p_Table)

!<description>
    ! This subroutine copies the content of the given integer array
    ! to the lists of the arraylist making use of the second table array
!</description>

!<input>
    ! pointer to the data
    integer, dimension(:), intent(in) :: p_IDataSrc

    ! pointer to the table
    integer, dimension(:), intent(in) :: p_Table
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: itable,ntable,ipos,kpos

    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListIntTbl')
      call sys_halt()
    end if
    
    ntable=size(p_Table)-1
    do itable=1,ntable
      do ipos=p_Table(itable),p_Table(itable+1)-1
        call arrlst_appendToArrayList(rarraylist,itable,p_IDataSrc(ipos),kpos)
      end do
    end do
  end subroutine arrlst_copyToArrayListIntTbl
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine arrlst_swapArrayList(rarraylist,itable,jtable)

!<description>
    ! This subroutine swaps two tables and the associated list in the
    ! given arraylist
!</description>
    
!<input>
    integer, intent(in) :: itable,jtable
!</input>

!<inputoutput>
    ! arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ihead,itail,iitem,ina
    
    ! Swap
    ihead = rarraylist%p_Ktable(ARRLST_HEAD,itable)
    itail = rarraylist%p_Ktable(ARRLST_TAIL,itable)
    iitem = rarraylist%p_Ktable(ARRLST_ITEM,itable)
    ina   = rarraylist%p_Ktable(ARRLST_NA,  itable)
    
    rarraylist%p_Ktable(ARRLST_HEAD,itable) = rarraylist%p_Ktable(ARRLST_HEAD,jtable)
    rarraylist%p_Ktable(ARRLST_TAIL,itable) = rarraylist%p_Ktable(ARRLST_TAIL,jtable)
    rarraylist%p_Ktable(ARRLST_ITEM,itable) = rarraylist%p_Ktable(ARRLST_ITEM,jtable)
    rarraylist%p_Ktable(ARRLST_NA,  itable) = rarraylist%p_Ktable(ARRLST_NA,  jtable)

    rarraylist%p_Ktable(ARRLST_HEAD,jtable) = ihead
    rarraylist%p_Ktable(ARRLST_TAIL,jtable) = itail
    rarraylist%p_Ktable(ARRLST_ITEM,jtable) = iitem
    rarraylist%p_Ktable(ARRLST_NA,  jtable) = ina
  end subroutine arrlst_swapArrayList

  ! ***************************************************************************
  
!<function>
  
  pure function arrlst_getFirstInArrayList(rarraylist,itable) result(ipos)

!<description>
    ! This function returns the position of the first list item in
    ! the given table. If the specified table does not exist than
    ! IPOS=-1 is returned.
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! Number of table
    integer, intent(in) :: itable
!</input>

!<result>
    ! position of first item
    integer :: ipos
!</result>
!</function>

    if (itable < 1 .or. itable > rarraylist%NTABLE) then
      ipos = -1
      return
    end if
    
    ipos=rarraylist%p_Knext(rarraylist%p_Ktable(ARRLST_HEAD,itable))
  end function arrlst_getFirstInArrayList

  ! ***************************************************************************
  
!<function>
  
  pure function arrlst_getLastInArrayList(rarraylist,itable) result(ipos)

!<description>
    ! This function returns the position of the last list item in the
    ! given table. Of the specified table does not exist than IPOS=-1
    ! is returned.
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! Number of table
    integer, intent(in) :: itable
!</input>

!<result>
    ! position of last item
    integer :: ipos
!</result>
!</function>

    if (itable < 1 .or. itable > rarraylist%NTABLE) then
      ipos = -1
      return
    end if
    
    ipos=rarraylist%p_Knext(rarraylist%p_Ktable(ARRLST_TAIL,itable))
  end function arrlst_getLastInArrayList

  ! ***************************************************************************

!<function>

  function arrlst_getNextInArrayList(rarraylist,itable,breset) result(ipos)

!<description>
    ! This function returns the position of the next free list item
    ! for a given table and resets the list if required. If the
    ! specified table does not exist then IPOS=-1 is returned.
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Reset list?
    logical, intent(in) :: breset
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! position of last item
    integer :: ipos
!</result>
!</function>
    
    if (itable < 1 .or. itable > rarraylist%NTABLE) then
      ipos = -1
      return
    end if

    ! Should we reset the item pointer?
    if (breset) then
      ipos = rarraylist%p_Ktable(ARRLST_HEAD,itable)
      rarraylist%p_Ktable(ARRLST_ITEM,itable) = ipos
      return
    end if

    ! Get next item and increase item pointer
    ipos = rarraylist%p_Knext(rarraylist%p_Ktable(ARRLST_ITEM,itable))
    rarraylist%p_Ktable(ARRLST_ITEM,itable) = ipos
  end function arrlst_getNextInArrayList

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_prependToArrayListDble(rarraylist,itable,da,ipos)

!<description>
    ! This subroutine prepends a Double data to the list of the given
    ! table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    real(DP), intent(in) :: da
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListDble')
      call sys_halt()
    end if
    
    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListDble')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)
    
    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
   
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_DData(ipos)                = da
    else
      rarraylist%p_Knext(ipos)                = rarraylist%p_Ktable(ARRLST_HEAD,itable)
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_DData(ipos)                = da
    end if
  end subroutine arrlst_prependToArrayListDble

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_prependToArrayListSngl(rarraylist,itable,sa,ipos)

!<description>
    ! This subroutine prepends a Single data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    real(SP), intent(in) :: sa
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListSngl')
      call sys_halt()
    end if
    
    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListSngl')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
     rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_FData(ipos)                = sa
    else
      rarraylist%p_Knext(ipos)                = rarraylist%p_Ktable(ARRLST_HEAD,itable)
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_FData(ipos)                = sa
    end if
  end subroutine arrlst_prependToArrayListSngl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_prependToArrayListInt(rarraylist,itable,ia,ipos)

!<description>
    ! This subroutine prepends an Integer data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    integer, intent(in) :: ia
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListInt')
      call sys_halt()
    end if
    
    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListInt')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
     rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_IData(ipos)                = ia
    else
      rarraylist%p_Knext(ipos)                = rarraylist%p_Ktable(ARRLST_HEAD,itable)
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_IData(ipos)                = ia
    end if
  end subroutine arrlst_prependToArrayListInt

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_appendToArrayListDble(rarraylist,itable,da,ipos)

!<description>
    ! This subroutine appends a Double data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    real(DP), intent(in) :: da
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListDble')
      call sys_halt()
    end if

    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListDble')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)
    
    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_DData(ipos)                = da
    else
      rarraylist%p_Knext(rarraylist%p_Ktable(ARRLST_TAIL,itable)) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable)                   = ipos
      rarraylist%p_Knext(ipos)                                  = ARRLST_NULL
      rarraylist%p_DData(ipos)                                  = da
    end if
  end subroutine arrlst_appendToArrayListDble

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_appendToArrayListSngl(rarraylist,itable,sa,ipos)

!<description>
    ! This subroutine appends a Single data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable
    
    ! Data
    real(SP), intent(in) :: sa
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListSngl')
      call sys_halt()
    end if

    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListSngl')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_FData(ipos)                = sa
    else
      rarraylist%p_Knext(rarraylist%p_Ktable(ARRLST_TAIL,itable)) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable)                   = ipos
      rarraylist%p_Knext(ipos)                                  = ARRLST_NULL
      rarraylist%p_FData(ipos)                                  = sa
    end if
  end subroutine arrlst_appendToArrayListSngl
  
  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_appendToArrayListInt(rarraylist,itable,ia,ipos)

!<description>
    ! This subroutine appends an Integer data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    integer, intent(in) :: ia
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the appended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListInt')
      call sys_halt()
    end if
    
    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListInt')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_IData(ipos)                = ia
    else
      rarraylist%p_Knext(rarraylist%p_Ktable(ARRLST_TAIL,itable)) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable)                   = ipos
      rarraylist%p_Knext(ipos)                                  = ARRLST_NULL
      rarraylist%p_IData(ipos)                                  = ia
    end if
  end subroutine arrlst_appendToArrayListInt

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_insertIntoArrayListDble(rarraylist,itable,da,ipred,ipos)

!<description>
    ! This subroutine inserts a new Double data into the list of a
    ! given table AFTER the position ipred
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    real(DP), intent(in) :: da

    ! Position of predecessor
    integer, intent(in) :: ipred
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListDble')
      call sys_halt()
    end if
    
    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListDble')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_DData(ipos)                = da
    elseif (ipred .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) then
      rarraylist%p_Knext(ipred)               = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_DData(ipos)                = da
    elseif (ipred < ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Knext(ipos)                = -ipred
      rarraylist%p_DData(ipos)                = da
    else
      rarraylist%p_Knext(ipos)                = rarraylist%p_Knext(ipred)
      rarraylist%p_Knext(ipred)               = ipos
      rarraylist%p_DData(ipos)                = da
    end if
  end subroutine arrlst_insertIntoArrayListDble

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_insertIntoArrayListSngl(rarraylist,itable,sa,ipred,ipos)

!<description>
    ! This subroutine inserts a new Single data into the list of a
    ! given table AFTER the position ipred
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    real(SP), intent(in) :: sa

    ! Position of predecessor
    integer, intent(in) :: ipred
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListSngl')
      call sys_halt()
    end if
    
    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListSngl')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos               = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_FData(ipos)                = sa
    elseif (ipred .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) then
      rarraylist%p_Knext(ipred)               = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_FData(ipos)                = sa
    elseif (ipred < ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Knext(ipos)                = -ipred
      rarraylist%p_FData(ipos)                = sa
    else
      rarraylist%p_Knext(ipos)                = rarraylist%p_Knext(ipred)
      rarraylist%p_Knext(ipred)               = ipos
      rarraylist%p_FData(ipos)                = sa
    end if
  end subroutine arrlst_insertIntoArrayListSngl

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_insertIntoArrayListInt(rarraylist,itable,ia,ipred,ipos)

!<description>
    ! This subroutine inserts a new Integer data into the list of a
    ! given table AFTER the position ipred
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    integer, intent(in) :: ia

    ! Position of predecessor
    integer, intent(in) :: ipred
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListInt')
      call sys_halt()
    end if
    
    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListInt')
      call sys_halt()
    end if
    if (rarraylist%NTABLE < itable) call arrlst_createArrayList(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ARRLST_NA,itable) =  rarraylist%p_Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%p_Knext(ARRLST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call arrlst_resizeArrayList(rarraylist,ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if
    
    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ARRLST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ARRLST_FREE) = rarraylist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    if (rarraylist%p_Ktable(ARRLST_HEAD,itable) .eq. ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_IData(ipos)                = ia
    elseif (ipred .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable)) then
      rarraylist%p_Knext(ipred)               = ipos
      rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)                = ARRLST_NULL
      rarraylist%p_IData(ipos)                = ia
    elseif (ipred < ARRLST_NULL) then
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Knext(ipos)                = -ipred
      rarraylist%p_IData(ipos)                = ia
    else
      rarraylist%p_Knext(ipos)                = rarraylist%p_Knext(ipred)
      rarraylist%p_Knext(ipred)               = ipos
      rarraylist%p_IData(ipos)                = ia
    end if
  end subroutine arrlst_insertIntoArrayListInt
  
  ! ***************************************************************************
  
!<function>

  function arrlst_deleteFromArrayListDble(rarraylist,itable,da) result(f)

!<description>
    ! This function deletes a Double data from the arraylist
    ! associated with the given table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    real(DP), intent(in) :: da
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / ARRAYLIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_deleteFromArrayListDble')
      call sys_halt()
    end if

    ! Search for data
    f=arrlst_searchInArrayList(rarraylist,itable,da,ipred)
    if (f .eq. ARRAYLIST_NOT_FOUND) return
    
    ! Delete data
    rarraylist%NA = rarraylist%NA-1
    rarraylist%p_Ktable(ARRLST_NA,itable) = rarraylist%p_Ktable(ARRLST_NA,itable)-1

    ! Are we first entry in list?
    if (ipred < 0) then

      ! Get position
      ipred = -ipred
      ipos  = rarraylist%p_Knext(ipred)
      
      ! Update free position
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Knext(ipred)               = rarraylist%p_Knext(ARRLST_FREE)
      rarraylist%p_Knext(ARRLST_FREE)         = -ipred

    else

      ! Get position
      ipos = rarraylist%p_Knext(ipred)
      if (rarraylist%p_Knext(ipred) .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable))&
          rarraylist%p_Ktable(ARRLST_TAIL,itable) = ipred
      
      ! Update free position
      rarraylist%p_Knext(ipred)       = rarraylist%p_Knext(ipos)
      rarraylist%p_Knext(ipos)        = rarraylist%p_Knext(ARRLST_FREE)
      rarraylist%p_Knext(ARRLST_FREE) = -ipos
    end if
  end function arrlst_deleteFromArrayListDble
  
  ! ***************************************************************************
  
!<function>

  function arrlst_deleteFromArrayListSngl(rarraylist,itable,sa) result(f)

!<description>
    ! This function deletes a Single data from the arraylist
    ! associated with the given table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    real(SP), intent(in) :: sa
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / ARRAYLIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_deleteFromArrayListSngl')
      call sys_halt()
    end if

    ! Search for data
    f=arrlst_searchInArrayList(rarraylist,itable,sa,ipred)
    if (f .eq. ARRAYLIST_NOT_FOUND) return

    ! Delete data
    rarraylist%NA = rarraylist%NA-1
    rarraylist%p_Ktable(ARRLST_NA,itable) = rarraylist%p_Ktable(ARRLST_NA,itable)-1
    
    ! Are we first entry in list?
    if (ipred < 0) then
      
      ! Get position
      ipred = -ipred
      ipos  = rarraylist%p_Knext(ipred)
      
      ! Update free position
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Knext(ipred)         = rarraylist%p_Knext(ARRLST_FREE)
      rarraylist%p_Knext(ARRLST_FREE)         = -ipred
      
    else
      
      ! Get position
      ipos = rarraylist%p_Knext(ipred)
      if (rarraylist%p_Knext(ipred) .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable))&
          rarraylist%p_Ktable(ARRLST_TAIL,itable)=ipred
      
      ! Update free position
      rarraylist%p_Knext(ipred) = rarraylist%p_Knext(ipos)
      rarraylist%p_Knext(ipos)  = rarraylist%p_Knext(ARRLST_FREE)
      rarraylist%p_Knext(ARRLST_FREE) = -ipos
    end if
  end function arrlst_deleteFromArrayListSngl

  ! ***************************************************************************
  
!<function>

  function arrlst_deleteFromArrayListInt(rarraylist,itable,ia) result(f)

!<description>
    ! This function deletes an Integer data from the arraylist
    ! associated with the given table.
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! Data
    integer, intent(in) :: ia
!</input>

!<inputoutput>
    ! list
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / ARRAYLIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_deleteFromArrayListInt')
      call sys_halt()
    end if

    ! Search for data
    f=arrlst_searchInArrayList(rarraylist,itable,ia,ipred)
    if (f .eq. ARRAYLIST_NOT_FOUND) return

    ! Delete data
    rarraylist%NA = rarraylist%NA-1
    rarraylist%p_Ktable(ARRLST_NA,itable) = rarraylist%p_Ktable(ARRLST_NA,itable)-1

    ! Are we first entry in list?
    if (ipred < 0) then

      ! Get position
      ipred = -ipred
      ipos  = rarraylist%p_Knext(ipred)

      ! Update free position
      rarraylist%p_Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%p_Knext(ipred)               = rarraylist%p_Knext(ARRLST_FREE)
      rarraylist%p_Knext(ARRLST_FREE)         = -ipred
      
    else

      ! Get position
      ipos = rarraylist%p_Knext(ipred)

      ! Check if last entry should be deleted
      if (rarraylist%p_Knext(ipred) .eq. rarraylist%p_Ktable(ARRLST_TAIL,itable))&
          rarraylist%p_Ktable(ARRLST_TAIL,itable)=ipred
      
      ! Update free position
      rarraylist%p_Knext(ipred) = rarraylist%p_Knext(ipos)
      rarraylist%p_Knext(ipos)  = rarraylist%p_Knext(ARRLST_FREE)
      rarraylist%p_Knext(ARRLST_FREE) = -ipos
    end if
  end function arrlst_deleteFromArrayListInt

  ! ***************************************************************************
  
!<function>

  function arrlst_searchInArrayListDble(rarraylist,itable,da,ipred) result(f)

!<description>
    ! This function searches for a given Double data in the list of a
    ! given table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! Data
    real(DP), intent(in) :: da
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer, intent(out) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOUT_FOUND / ARRAYLIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ihead,itail,inext

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_searchInArrayListDble')
      call sys_halt()
    end if

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    if (itable < 1 .or. itable > rarraylist%NTABLE) return
    ipred=-rarraylist%p_Ktable(ARRLST_HEAD,itable)
    
    ! Check if list is empty
    if (ipred .eq. ARRLST_NULL) return

    ! Initialization
    ihead=rarraylist%p_Ktable(ARRLST_HEAD,itable)
    itail=rarraylist%p_Ktable(ARRLST_TAIL,itable)

    ! What kind of ordering are we
    select case(rarraylist%cordering)
    case (ARRAYLIST_UNORDERED)

      ! Check first item separately
      if (rarraylist%p_DData(ihead) .eq. da) then
        f=ARRAYLIST_FOUND; return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_DData(inext) .eq. da) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (inext .eq. itail) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_INCREASING)

      ! Check first item separately
      if (rarraylist%p_DData(ihead) .eq. da) then
        f=ARRAYLIST_FOUND; return
      elseif(rarraylist%p_DData(ihead) > da) then
        return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_DData(inext) .eq. da) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_DData(inext) > da) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_DECREASING)

      ! Check first item separately
      if (rarraylist%p_DData(ihead) .eq. da) then
        f=ARRAYLIST_FOUND; return
      elseif(rarraylist%p_DData(ihead) < da) then
        return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_DData(inext) .eq. da) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_DData(inext) < da) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_CSR7)

      ! Check first item separately
      if (rarraylist%p_DData(ihead) .eq. da) then
        f=ARRAYLIST_FOUND; return
      elseif(abs(itable-da) .le. SYS_EPSREAL_DP) then
        return
      end if
      ipred=-ipred
      
      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_DData(inext) .eq. da) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_DData(inext) > da) exit
        if (rarraylist%p_Knext(ipred) .eq. itail) then
          ipred=rarraylist%p_Knext(ipred); exit
        end if
        ipred=rarraylist%p_Knext(ipred)
      end do
    end select
  end function arrlst_searchInArrayListDble
  
  ! ***************************************************************************
  
!<function>

  function arrlst_searchInArrayListSngl(rarraylist,itable,sa,ipred) result(f)

!<description>
    ! This function searches for a given Single data in the list of a
    ! given table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! Data
    real(SP), intent(in) :: sa
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer, intent(out) :: ipred
!</output>

!<result>
    ! Result of the search ARRAYLIST_NOT_FOUND / ARRAYLIST_FOUND
    integer :: f
!</result>
!</function>
    
    ! local variables
    integer :: ihead,itail,inext

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_searchInArrayListSngl')
      call sys_halt()
    end if

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    if (itable < 1 .or. itable > rarraylist%NTABLE) return
    ipred=-rarraylist%p_Ktable(ARRLST_HEAD,itable)

    ! Check if list is empty
    if (ipred .eq. ARRLST_NULL) return

    ! Initialization
    ihead=rarraylist%p_Ktable(ARRLST_HEAD,itable)
    itail=rarraylist%p_Ktable(ARRLST_TAIL,itable)

    ! What kind of ordering are we
    select case(rarraylist%cordering)
    case (ARRAYLIST_UNORDERED)

      ! Check first item separately
      if (rarraylist%p_FData(ihead) .eq. sa) then
        f=ARRAYLIST_FOUND; return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_FData(inext) .eq. sa) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (inext .eq. itail) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_INCREASING)

      ! Check first item separately
      if (rarraylist%p_FData(ihead) .eq. sa) then
        f=ARRAYLIST_FOUND; return
      elseif(rarraylist%p_FData(ihead) > sa) then
        return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_FData(inext) .eq. sa) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_FData(inext) > sa) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_DECREASING)

      ! Check first item separately
      if (rarraylist%p_FData(ihead) .eq. sa) then
        f=ARRAYLIST_FOUND; return
      elseif(rarraylist%p_FData(ihead) < sa) then
        return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_FData(inext) .eq. sa) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_FData(inext) < sa) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_CSR7)

      ! Check first item separately
      if (rarraylist%p_FData(ihead) .eq. sa) then
        f=ARRAYLIST_FOUND; return
      elseif(abs(itable-sa) .le. SYS_EPSREAL_DP) then
        return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_FData(inext) .eq. sa) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_FData(inext) > sa) exit
        if (rarraylist%p_Knext(ipred) .eq. itail) then
          ipred=rarraylist%p_Knext(ipred); exit
        end if
        ipred=rarraylist%p_Knext(ipred)
      end do
    end select
  end function arrlst_searchInArrayListSngl

  ! ***************************************************************************
  
!<function>

  function arrlst_searchInArrayListInt(rarraylist,itable,ia,ipred) result(f)

!<description>
    ! This function searches for a given Integer data in the list of
    ! a given table
!</description>

!<input>
    ! Number of table
    integer, intent(in) :: itable

    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! Data
    integer, intent(in) :: ia
!</input>

!<output>
    ! Position of the predecessor of the found item
    integer, intent(out) :: ipred
!</output>

!<result>
    ! Result of the search ARRAYLIST_NOT_FOUND / ARRAYLIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ihead,itail,inext

    ! Check if list format is ok
    if (rarraylist%carraylistFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_searchInArrayListInt')
      call sys_halt()
    end if

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    if (itable < 1 .or. itable > rarraylist%NTABLE) return
    ipred=-rarraylist%p_Ktable(ARRLST_HEAD,itable)
    
    ! Check if list is empty
    if (ipred .eq. ARRLST_NULL) return

    ! Initialization
    ihead=rarraylist%p_Ktable(ARRLST_HEAD,itable)
    itail=rarraylist%p_Ktable(ARRLST_TAIL,itable)
        
    ! What kind of ordering are we
    select case(rarraylist%cordering)
    case (ARRAYLIST_UNORDERED)

      ! Check first item separately
      if (rarraylist%p_IData(ihead) .eq. ia) then
        f=ARRAYLIST_FOUND; return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_IData(inext) .eq. ia) then
          f=ARRAYLIST_FOUND; exit
        end if

        if (inext .eq. itail) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_INCREASING)

      ! Check first item separately
      if (rarraylist%p_IData(ihead) .eq. ia) then
        f=ARRAYLIST_FOUND; return
      elseif(rarraylist%p_IData(ihead) > ia) then
        return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_IData(inext) .eq. ia) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_IData(inext) > ia) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_DECREASING)

      ! Check first item separately
      if (rarraylist%p_IData(ihead) .eq. ia) then
        f=ARRAYLIST_FOUND; return
      elseif(rarraylist%p_IData(ihead) < ia) then
        return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_IData(inext) .eq. ia) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_IData(inext) < ia) exit
        ipred=rarraylist%p_Knext(ipred)
      end do
      
    case (ARRAYLIST_CSR7)

      ! Check first item separately
      if (rarraylist%p_IData(ihead) .eq. ia) then
        f=ARRAYLIST_FOUND; return
      elseif(itable .eq. ia) then
        return
      end if
      ipred=-ipred

      do while(ipred.ne.itail)
        inext = rarraylist%p_Knext(ipred)
        if (rarraylist%p_IData(inext) .eq. ia) then
          f=ARRAYLIST_FOUND; exit
        end if
        
        if (rarraylist%p_IData(inext) > ia) exit
        if (rarraylist%p_Knext(ipred) .eq. itail) then
          ipred=rarraylist%p_Knext(ipred); exit
        end if
        ipred=rarraylist%p_Knext(ipred)
      end do
    end select
  end function arrlst_searchInArrayListInt

  ! ***************************************************************************
  
!<subroutine>

  subroutine arrlst_printArrayList(rarraylist,itable)

!<description>
    ! This subroutine prints the content of the list
!</description>

!<input>
    ! list
    type(t_arraylist), intent(in) :: rarraylist

    ! OPTIONAL: number of table
    integer, intent(in), optional :: itable
!</input>
!</subroutine>

    ! local variable
    integer :: iitable,itable1,itable2,ipos,itail
    
    if (present(itable)) then
      itable1=itable; itable2=itable
    else
      itable1=1; itable2=rarraylist%NTABLE
    end if

    select case (rarraylist%carraylistFormat)
    case (ST_DOUBLE)

      ! Loop over all selected tables
      do iitable=itable1,itable2
        call output_line('Table number: '//trim(sys_siL(iitable,15)))
        call output_line('----------------------------------------')

        ipos  = rarraylist%p_Ktable(ARRLST_HEAD,iitable)
        if (ipos .eq. ARRLST_NULL) cycle
        itail = rarraylist%p_Ktable(ARRLST_TAIL,iitable)
        
        do
          write(*,FMT='(A,",")',ADVANCE='NO') trim(sys_sdL(rarraylist%p_DData(ipos),8))
          if (ipos .eq. itail) exit
          ipos = rarraylist%p_Knext(ipos)
        end do
      end do
      
    case (ST_SINGLE)

      ! Loop over all selected tables
      do iitable=itable1,itable2
        call output_line('Table number: '//trim(sys_siL(iitable,15)))
        call output_line('----------------------------------------')
        
        ipos = rarraylist%p_Ktable(ARRLST_HEAD,iitable)
        if (ipos .eq. ARRLST_NULL) cycle
        itail = rarraylist%p_Ktable(ARRLST_TAIL,iitable)
        
        do
          write(*,FMT='(A,",")',ADVANCE='NO') trim(sys_sdL(real(rarraylist%p_FData(ipos),DP),8))
          if (ipos .eq. itail) exit
          ipos = rarraylist%p_Knext(ipos)
        end do
      end do
      
    case (ST_INT)
      
      ! Loop over all selected tables
      do iitable=itable1,itable2
        call output_line('Table number: '//trim(sys_siL(iitable,15)))
        call output_line('----------------------------------------')

        ipos = rarraylist%p_Ktable(ARRLST_HEAD,iitable)
        if (ipos .eq. ARRLST_NULL) cycle
        itail = rarraylist%p_Ktable(ARRLST_TAIL,iitable)
        
        do
          write(*,FMT='(A,",")',ADVANCE='NO') trim(sys_siL(rarraylist%p_IData(ipos),15))
          if (ipos .eq. itail) exit
          ipos = rarraylist%p_Knext(ipos)
        end do
        write(*,*)
      end do
      
    case DEFAULT
      call output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_printArrayList')
      call sys_halt()
    end select
  end subroutine arrlst_printArrayList

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_infoArrayList(rarraylist)

!<description>
    ! This subroutine prints information about the arraylist
!</description>

!<input>
    ! arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>
!</subroutine>

    call output_line('Arraylist:')
    call output_line('----------')
    call output_line('NA:       '//trim(sys_siL(rarraylist%NA,15)))
    call output_line('NNA:      '//trim(sys_siL(rarraylist%NNA,15)))
    call output_line('NNA0:     '//trim(sys_siL(rarraylist%NNA0,15)))
    call output_line('NTABLE:   '//trim(sys_siL(rarraylist%NTABLE,15)))
    call output_line('NNTABLE:  '//trim(sys_siL(rarraylist%NNTABLE,15)))
    call output_line('NNTABLE0: '//trim(sys_siL(rarraylist%NNTABLE0,15)))
    call output_line('NRESIZE:  '//trim(sys_siL(rarraylist%NRESIZE,15)))
    call output_line('dfactor:  '//trim(sys_sdL(rarraylist%dfactor,2)))
    call output_line('h_Ktable: '//trim(sys_siL(rarraylist%h_Ktable,15)))
    call output_line('h_Knext:  '//trim(sys_siL(rarraylist%h_Knext,15)))
    call output_line('h_Data:   '//trim(sys_siL(rarraylist%h_Data,15)))
    call output_lbrk()
    call output_line('Current data  memory usage: '//&
        trim(sys_sdL(100*rarraylist%NA/&
        real(max(1,rarraylist%NNA),DP),2))//'%')
    call output_line('Current table memory usage: '//&
        trim(sys_sdL(100*rarraylist%NTABLE/&
        real(max(1,rarraylist%NNTABLE),DP),2))//'%')
    call output_line('Total data memory usage:    '//&
        trim(sys_sdL(100*rarraylist%NNA/&
        real(max(1,rarraylist%NNA0),DP),2))//'%')
    call output_line('Total table memory usage:   '//&
        trim(sys_sdL(100*rarraylist%NNTABLE/&
        real(max(1,rarraylist%NNTABLE0),DP),2))//'%')
    call output_lbrk()
  end subroutine arrlst_infoArrayList

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_duplicateArrayList(rarraylist,rarraylistBackup)

!<description>
    ! This subroutine makes a copy of an arraylist in memory.
    ! It does not make sense to share some information between arraylists,
    ! so each vectors is physically copied from the source arraylist
    ! to the destination arraylist.
!</description>

!<input>
    ! Source arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<inputoutput>
    ! Destination arraylist
    type(t_arraylist), intent(inout) :: rarraylistBackup
!</inputoutput>
!</subroutine>

    ! Release backup arraylist
    call arrlst_releaseArrayList(rarraylistBackup)

    ! Copy all data
    rarraylistBackup = rarraylist

    ! Reset handles
    rarraylistBackup%h_Ktable = ST_NOHANDLE
    rarraylistBackup%h_Knext  = ST_NOHANDLE
    rarraylistBackup%h_Data   = ST_NOHANDLE

    ! Copy storage blocks
    if (rarraylist%h_Ktable .ne. ST_NOHANDLE) then
      call storage_copy(rarraylist%h_Ktable,&
          rarraylistBackup%h_Ktable)
      call storage_getbase_int2D(rarraylistBackup%h_Ktable,&
          rarraylistBackup%p_Ktable)
    end if

    if (rarraylist%h_Knext .ne. ST_NOHANDLE) then
      call storage_copy(rarraylist%h_Knext,&
          rarraylistBackup%h_Knext)
      call storage_getbase_int(rarraylistBackup%h_Knext,&
          rarraylistBackup%p_Knext)
    end if

    if (rarraylist%h_Data .ne. ST_NOHANDLE) then
      call storage_copy(rarraylist%h_Data,&
          rarraylistBackup%h_Data)

      select case(rarraylistBackup%carraylistFormat)
      case(ST_DOUBLE)
        call storage_getbase_double(rarraylistBackup%h_Data,&
            rarraylistBackup%p_DData)
      case(ST_SINGLE)
        call storage_getbase_single(rarraylistBackup%h_Data,&
            rarraylistBackup%p_FData)
      case(ST_INT)
        call storage_getbase_int(rarraylistBackup%h_Data,&
            rarraylistBackup%p_IData)
      case DEFAULT
        call output_line('Unsupported data format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'arrlst_duplicateArrayList')
        call sys_halt()
      end select
    end if
  end subroutine arrlst_duplicateArrayList

  ! ***************************************************************************

!<subroutine>

  subroutine arrlst_restoreArrayList(rarraylistBackup,rarraylist)

!<description>
    ! This subroutine restores an arraylist from a previous backup.
    ! The format and ordering of both arraylists must be the same.
!</description>

!<input>
    ! Backup of an arraylist
    type(t_arraylist), intent(in) :: rarraylistBackup
!</input>

!<inputoutput>
    ! Destination arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Check that both arraylists are compatible
    if (rarraylist%carraylistFormat .ne. rarraylistBackup%carraylistFormat .or.&
        rarraylist%cordering        .ne. rarraylistBackup%cordering) then
      call output_line('Incompatible arraylists!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_restoreArrayList')
      call sys_halt()
    end if

    ! Release arraylist
    call arrlst_releaseArrayList(rarraylist)

    ! Duplicate the backup
    call arrlst_duplicateArrayList(rarraylistBackup,rarraylist)
  end subroutine arrlst_restoreArrayList
end module arraylist
