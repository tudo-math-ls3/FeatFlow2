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
!# 1.) arrlst_createArrayList
!#     -> Create an empty arraylist
!#
!# 2.) arrlst_releaseArrayList = arrlst_releaseArrayList /
!#                               arrlst_releaseArrayList_table
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
!# 5.) arrlst_copyArrayListTable = arrlst_copyFromArrayList_table /
!#                                 arrlst_copyFromArrayListDble_table /
!#                                 arrlst_copyFromArrayListSngl_table /
!#                                 arrlst_copyFromArrayListInt_table /
!#                                 arrlst_copyToArrayList_table /
!#                                 arrlst_copyToArrayListDble_table /
!#                                 arrlst_copyToArrayListSngl_table /
!#                                 arrlst_copyToArrayListInt_table
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
!# </purpose>
!##############################################################################

MODULE arraylist
  USE fsystem
  USE storage
  USE genoutput
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_arraylist
  PUBLIC :: arrlst_createArrayList
  PUBLIC :: arrlst_releaseArrayList
  PUBLIC :: arrlst_resizeArrayList
  PUBLIC :: arrlst_copyArrayList
  PUBLIC :: arrlst_copyArrayListTable
  PUBLIC :: arrlst_swapArrayList
  PUBLIC :: arrlst_getFirstInArrayList
  PUBLIC :: arrlst_getLastInArrayList
  PUBLIC :: arrlst_getNextInArrayList
  PUBLIC :: arrlst_prependToArrayList
  PUBLIC :: arrlst_appendToArrayList
  PUBLIC :: arrlst_insertIntoArrayList
  PUBLIC :: arrlst_deleteFromArrayList
  PUBLIC :: arrlst_searchInArrayList
  PUBLIC :: arrlst_printArrayList
  PUBLIC :: arrlst_infoArrayList

!<constants>

!<constantblock description="KIND values for list data">

  ! kind value for indices in arraylist
  INTEGER, PARAMETER, PUBLIC :: PREC_ARRAYLISTIDX = I32

  ! kind value for indices in table
  INTEGER, PARAMETER, PUBLIC :: PREC_TABLEIDX = I32

!</constantblock>

!<constantblock description="Global flags for arraylist ordering">

  ! Identifier for unordered arraylist
  INTEGER, PARAMETER, PUBLIC :: ARRAYLIST_UNORDERED  = 0
  
  ! Identifier for increasingly ordered list
  INTEGER, PARAMETER, PUBLIC :: ARRAYLIST_INCREASING = 1

  ! Identifier for decreasingly ordered arraylist
  INTEGER, PARAMETER, PUBLIC :: ARRAYLIST_DECREASING = 2

  ! Identifier for ordered arraylist 
  INTEGER, PARAMETER, PUBLIC :: ARRAYLIST_CSR7       = 3

!</constantblock>

!<constantblock description="Global flags for arraylist operations">

  ! Identifier for "not found in arraylist"
  INTEGER, PARAMETER, PUBLIC :: ARRAYLIST_NOT_FOUND = -1

  ! Identifier for "found in arraylist"
  INTEGER, PARAMETER, PUBLIC :: ARRAYLIST_FOUND     =  0
  
!</constantblock>

!<constantblock description="Internal tags for arraylist status">
  
  ! Tag for empty arraylist
  INTEGER, PARAMETER, PUBLIC :: ARRLST_NULL =  0

  ! Tag for next free position in storage of arraylist
  INTEGER, PARAMETER :: ARRLST_FREE = 0

  ! Tag for head of each list
  INTEGER, PARAMETER :: ARRLST_HEAD = 1

  ! Tag for tail of each list
  INTEGER, PARAMETER :: ARRLST_TAIL = 2

  ! Tag for last item stored in each list
  INTEGER, PARAMETER :: ARRLST_ITEM = 3

  ! Tag for number of entries stored in each list
  INTEGER, PARAMETER :: ARRLST_NA   = 4

!</constantblock>
!</constants>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<types>
!<typeblock>

  TYPE t_arraylist
    ! Format-Tag:
    INTEGER :: carraylistFormat                   = ST_NOHANDLE
    
    ! Type of arraylist ordering
    INTEGER :: cordering                          = ARRAYLIST_UNORDERED
    
    ! Number of tables that are currently stored in the arraylist
    INTEGER(PREC_TABLEIDX) :: NTABLE              = 0
    
    ! Total number of tables that can be stored in the arraylist
    INTEGER(PREC_TABLEIDX) :: NNTABLE             = 0

    ! Total number of tables that can initially be stored in the
    ! arraylist. This information is needed to compute the growth of
    ! the arraylist after several resize operations
    INTEGER(PREC_TABLEIDX) :: NNTABLE0            = 0

    ! Number of items that are currently stored in all lists
    INTEGER(PREC_ARRAYLISTIDX) :: NA              = 0

    ! Total number of items that can be stored in all lists
    INTEGER(PREC_ARRAYLISTIDX) :: NNA             = 0

    ! Total number of items that can initially be stored in all
    ! lists. This information is needed to compute the growth of the
    ! arraylist after several resize operations
    INTEGER(PREC_ARRAYLISTIDX) :: NNA0            = 0

    ! Number of resize operations performed with the arraylist
    INTEGER :: NRESIZE                            = 0

    ! Factor by which the arraylist is enlarged if new storage is allocate
    REAL(DP) :: dfactor                           = 1.5_DP

    ! Handle to the lookup table
    INTEGER :: h_Ktable                           = ST_NOHANDLE

    ! Handle to the arraylist structure
    INTEGER :: h_Knext                            = ST_NOHANDLE
    
    ! Handle to the arraylist data
    INTEGER :: h_Data                             = ST_NOHANDLE

    ! Table structure
    ! NOTE: This array is introduced to increase performance. It
    ! should not be touched by the user. However, if the handle would
    ! be dereferenced for each operation such as search, delete,
    ! performance would be very poor.
    INTEGER(PREC_ARRAYLISTIDX), DIMENSION(:,:), POINTER :: Ktable => NULL()
    
    ! ArrayList structure
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_ARRAYLISTIDX), DIMENSION(:), POINTER :: Knext => NULL()

    ! ArrayList data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(DP), DIMENSION(:), POINTER :: DData => NULL()

    ! ArrayList data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(SP), DIMENSION(:), POINTER :: SData => NULL()

    ! ArrayList data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_ARRAYLISTIDX), DIMENSION(:), POINTER :: IData => NULL()
  END TYPE t_arraylist
  
!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************
  
  INTERFACE arrlst_releaseArrayList
    MODULE PROCEDURE arrlst_releaseArrayList
    MODULE PROCEDURE arrlst_releaseArrayList_table
  END INTERFACE 
  
  INTERFACE arrlst_copyArrayList
    MODULE PROCEDURE arrlst_copyFromArrayList
    MODULE PROCEDURE arrlst_copyFromArrayListDble
    MODULE PROCEDURE arrlst_copyFromArrayListSngl
    MODULE PROCEDURE arrlst_copyFromArrayListInt
    MODULE PROCEDURE arrlst_copyToArrayList
    MODULE PROCEDURE arrlst_copyToArrayListDble
    MODULE PROCEDURE arrlst_copyToArrayListSngl
    MODULE PROCEDURE arrlst_copyToArrayListInt
  END INTERFACE

  INTERFACE arrlst_copyArrayListTable
    MODULE PROCEDURE arrlst_copyFromArrayList_table
    MODULE PROCEDURE arrlst_copyFromArrayListDble_table
    MODULE PROCEDURE arrlst_copyFromArrayListSngl_table
    MODULE PROCEDURE arrlst_copyFromArrayListInt_table
    MODULE PROCEDURE arrlst_copyToArrayList_table
    MODULE PROCEDURE arrlst_copyToArrayListDble_table
    MODULE PROCEDURE arrlst_copyToArrayListSngl_table
    MODULE PROCEDURE arrlst_copyToArrayListInt_table
  END INTERFACE
    
  INTERFACE arrlst_prependToArrayList
    MODULE PROCEDURE arrlst_prependToArrayListDble
    MODULE PROCEDURE arrlst_prependToArrayListSngl
    MODULE PROCEDURE arrlst_prependToArrayListInt
  END INTERFACE
   
  INTERFACE arrlst_appendToArrayList
    MODULE PROCEDURE arrlst_appendToArrayListDble
    MODULE PROCEDURE arrlst_appendToArrayListSngl
    MODULE PROCEDURE arrlst_appendToArrayListInt
  END INTERFACE
   
  INTERFACE arrlst_insertIntoArrayList
    MODULE PROCEDURE arrlst_insertIntoArrayListDble
    MODULE PROCEDURE arrlst_insertIntoArrayListSngl
    MODULE PROCEDURE arrlst_insertIntoArrayListInt
  END INTERFACE
 
  INTERFACE arrlst_deleteFromArrayList
    MODULE PROCEDURE arrlst_deleteFromArrayListDble
    MODULE PROCEDURE arrlst_deleteFromArrayListSngl
    MODULE PROCEDURE arrlst_deleteFromArrayListInt
  END INTERFACE
    
  INTERFACE arrlst_searchInArrayList
    MODULE PROCEDURE arrlst_searchInArrayListDble
    MODULE PROCEDURE arrlst_searchInArrayListSngl
    MODULE PROCEDURE arrlst_searchInArrayListInt
  END INTERFACE
  
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

CONTAINS
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE arrlst_createArrayList(rarraylist,nntable,nna,carraylistFormat,cordering,dfactor)

!<description>
    ! This subroutine creates a new arraylist
!</description>

!<input>
    ! Total number of tables
    INTEGER(PREC_TABLEIDX), INTENT(IN)     :: nntable

    ! Total number of items that can be stored in all lists
    INTEGER(PREC_ARRAYLISTIDX), INTENT(IN) :: nna

    ! Format-tag. Type of list format (Double,Single,Integer)
    INTEGER, INTENT(IN)                    :: carraylistFormat

    ! OPTIONAL: Format-tag. Type of list ordering 
    INTEGER, INTENT(IN), OPTIONAL          :: cordering
    
    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    REAL(DP), INTENT(IN), OPTIONAL         :: dfactor
!</input>

!<output>
    ! list
    TYPE(t_arraylist), INTENT(OUT)         :: rarraylist
!</output>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    
    ! Set factor
    IF (PRESENT(dfactor)) THEN
      IF (dfactor > 1_DP) rarraylist%dfactor=dfactor
    END IF

    ! Set list format
    rarraylist%carraylistFormat=carraylistFormat

    ! Set ordering
    IF (PRESENT(cordering)) THEN
      rarraylist%cordering=cordering
    END IF

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
    CALL storage_new('arrlst_createArrayList','Ktable',Isize,&
        ST_INT,rarraylist%h_Ktable,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int2D(rarraylist%h_Ktable,rarraylist%Ktable)
    
    CALL storage_new('arrlst_createArrayList','Knext',ARRLST_FREE,nna,ST_INT,&
        rarraylist%h_Knext,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(rarraylist%h_Knext,rarraylist%Knext)

    SELECT CASE(rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_new('arrlst_createArrayList','Data',nna,ST_DOUBLE,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double(rarraylist%h_Data,rarraylist%DData)
      
    CASE (ST_SINGLE)
      CALL storage_new('arrlst_createArrayList','Data',nna,ST_SINGLE,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_single(rarraylist%h_Data,rarraylist%SData)
      
    CASE (ST_INT)
      CALL storage_new('arrlst_createArrayList','Data',nna,ST_INT,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int(rarraylist%h_Data,rarraylist%IData)
      
    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_createArrayList')
      CALL sys_halt()
    END SELECT
    
    ! Initialize list structures
    rarraylist%Knext(ARRLST_FREE)    = 1
  END SUBROUTINE arrlst_createArrayList
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_createArrayList_table(rarraylist,itable)

!<description>
    ! This subroutine creates new table entries up to position
    ! itable. Note: This subroutine cannot be called from outside of
    ! this module by the user. It is only used internally if some
    ! item should be added to a list for which there exists no table.
!</description>

!<input>
    ! Number of last table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! Check if table number is valid
    IF (itable < 1) THEN
      CALL output_line('Invalid table number',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_createArrayList_table')
      CALL sys_halt()
    END IF

    ! Resize tables if required
    IF (itable > rarraylist%NNTABLE) CALL arrlst_resizeArrayList_table(&
        rarraylist,CEILING(itable*rarraylist%dfactor))
    
    ! Initialize structures
    rarraylist%Ktable(ARRLST_HEAD:ARRLST_NA,rarraylist%NTABLE+1:itable) = ARRLST_NULL

    ! Set new table size
    rarraylist%NTABLE = MAX(rarraylist%NTABLE,itable)
  END SUBROUTINE arrlst_createArrayList_table

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE arrlst_releaseArrayList(rarraylist)

!<description>
    ! This subroutine releases an existing arraylist
!</description>

!<inputoutput>
    TYPE(t_arraylist), INTENT(INOUT) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Release memory
    IF (rarraylist%h_Ktable .NE. ST_NOHANDLE) CALL storage_free(rarraylist%h_Ktable)
    IF (rarraylist%h_Knext .NE. ST_NOHANDLE)  CALL storage_free(rarraylist%h_Knext)
    IF (rarraylist%h_Data .NE. ST_NOHANDLE)   CALL storage_free(rarraylist%h_Data)
    NULLIFY(rarraylist%Ktable,rarraylist%Knext,rarraylist%Ddata,&
        rarraylist%SData,rarraylist%IData)

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
  END SUBROUTINE arrlst_releaseArrayList

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_releaseArrayList_table(rarraylist,itable)

!<description>
    ! This subroutine releases a table from the arraylist
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<inputoutput>
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>
!</subroutine>

    ! Check if table exists
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_releaseArrayList_table')
      CALL sys_halt()
    END IF

    ! Decrease number of entries by the number of entries present 
    ! in the table which is released
    rarraylist%NA = rarraylist%NA - rarraylist%Ktable(ARRLST_NA,itable)

    ! Reset table
    rarraylist%Ktable(ARRLST_HEAD:ARRLST_NA,itable) = ARRLST_NULL

    ! Decrease number of tables if the last table has been deleted
    IF (itable .EQ. rarraylist%NTABLE)&
        rarraylist%NTABLE = rarraylist%NTABLE-1
  END SUBROUTINE arrlst_releaseArrayList_table

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE arrlst_resizeArrayList(rarraylist,nna)

!<description>
    ! This subroutine reallocates memory for an existing list
!</description>

!<input>
    ! New number of total items that can be stored in the list
    INTEGER(PREC_ARRAYLISTIDX), INTENT(IN) :: nna
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)       :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! Set new size and increase counter
    rarraylist%NNA     = nna
    rarraylist%NRESIZE = rarraylist%NRESIZE+1

    CALL storage_realloc('arrlst_resizeArrayList',ARRLST_FREE,nna,&
        rarraylist%h_Knext,ST_NEWBLOCK_NOINIT,.TRUE.)
    CALL storage_realloc('arrlst_resizeArrayList',nna,&
        rarraylist%h_Data,ST_NEWBLOCK_NOINIT,.TRUE.)
    CALL storage_getbase_int(rarraylist%h_Knext,rarraylist%Knext)

    SELECT CASE(rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(rarraylist%h_Data,rarraylist%DData)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(rarraylist%h_Data,rarraylist%SData)

    CASE (ST_INT)
      CALL storage_getbase_int(rarraylist%h_Data,rarraylist%IData)

    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_resizeArrayList')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE arrlst_resizeArrayList

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_resizeArrayList_table(rarraylist,nntable)

!<description>
    ! This subroutine reallocates memory for the lookup table
!</description>

!<input>
    ! New number of tables
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: nntable
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>
!</subroutine>

    ! Set new size
    rarraylist%NNTABLE = nntable
    rarraylist%NRESIZE = rarraylist%NRESIZE+1

    CALL storage_realloc('arrlst_resizeArrayList_table',nntable,&
        rarraylist%h_Ktable,ST_NEWBLOCK_NOINIT,.TRUE.)
    CALL storage_getbase_int2D(rarraylist%h_Ktable,rarraylist%Ktable)

!!$    It should not be necessary to clear all arrays
!!$    rarraylist%Ktable(ARRLST_HEAD,rarraylist%NTABLE+1:) = ARRLST_NULL
!!$    rarraylist%Ktable(ARRLST_TAIL,rarraylist%NTABLE+1:) = ARRLST_NULL
!!$    rarraylist%Ktable(ARRLST_ITEM,rarraylist%NTABLE+1:) = ARRLST_HEAD
!!$    rarraylist%Ktable(ARRLST_NA,  rarraylist%NTABLE+1:) = ARRLST_NULL
  END SUBROUTINE arrlst_resizeArrayList_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyFromArrayList(rarraylist,itable,h_Data)

!<description>
    ! This subroutine copies the content of the list of a given table
    ! to the given handle.
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)      :: rarraylist

    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<inputoutput>
    ! handle to the data
    INTEGER, INTENT(INOUT)             :: h_Data
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_DData
    REAL(SP), DIMENSION(:), POINTER :: p_SData
    INTEGER,  DIMENSION(:), POINTER :: p_IData
    INTEGER(I32) :: isize

    ! Transform the content of the list to h_Data
    IF (h_Data .EQ. ST_NOHANDLE) THEN
      CALL storage_new('arrlst_copyFromArrayList','Data',rarraylist%NA,&
          rarraylist%carraylistFormat,h_Data,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Data,isize)
      IF (isize < rarraylist%NA) THEN
        CALL storage_realloc('arrlst_copyFromArrayList',rarraylist%NA,&
            h_Data,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF
    
    ! What kind of data are we?
    SELECT CASE(rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_Data,p_DData)
      CALL arrlst_copyArrayList(rarraylist,itable,p_DData)
      
    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Data,p_SData)
      CALL arrlst_copyArrayList(rarraylist,itable,p_SData)
      
    CASE (ST_INT)
      CALL storage_getbase_int(h_Data,p_IData)
      CALL arrlst_copyArrayList(rarraylist,itable,p_IData)
      
    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayList')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE arrlst_copyFromArrayList

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyFromArrayList_table(rarraylist,h_Data,h_Table)

!<description>
    ! This subroutine copies the content of the table and all lists
    ! to the given handles.
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN) :: rarraylist

!</input>

!<inputoutput>
    ! handle to the data
    INTEGER, INTENT(INOUT)        :: h_Data

    ! handle to the table
    INTEGER, INTENT(INOUT)        :: h_Table
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TABLEIDX), DIMENSION(:), POINTER :: p_Table
    REAL(DP), DIMENSION(:), POINTER :: p_DData
    REAL(SP), DIMENSION(:), POINTER :: p_SData
    INTEGER,  DIMENSION(:), POINTER :: p_IData
    INTEGER(I32) :: isize

    ! Transform the content of the arraylist and the 
    ! table to h_Data and h_Table, respectively
    IF (h_Table .EQ. ST_NOHANDLE) THEN
      CALL storage_new('arrlst_copyFromArrayList_table','Table',&
          rarraylist%NTABLE+1,ST_INT,h_Table,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Table,isize)
      IF (isize < rarraylist%NTABLE+1) THEN
        CALL storage_realloc('arrlst_copyFromArrayList_table',&
            rarraylist%NTABLE+1,h_Table,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF
    CALL storage_getbase_int(h_Table,p_Table)
    
    IF (h_Data .EQ. ST_NOHANDLE) THEN
      CALL storage_new('arrlst_copyFromArrayList_table','Data',&
          rarraylist%NA,rarraylist%carraylistFormat,h_Data,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Data,isize)
      IF (isize < rarraylist%NA) THEN
        CALL storage_realloc('arrlst_copyFromArrayList_table',&
            rarraylist%NA,h_Data,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF

    ! What kind of array are we?
    SELECT CASE(rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_Data,p_DData)
      CALL arrlst_copyArrayListTable(rarraylist,p_DData,p_Table)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Data,p_SData)
      CALL arrlst_copyArrayListTable(rarraylist,p_SData,p_Table)

    CASE (ST_INT)
      CALL storage_getbase_int(h_Data,p_IData)
      CALL arrlst_copyArrayListTable(rarraylist,p_IData,p_Table)
      
    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayList_table')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE arrlst_copyFromArrayList_table
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyFromArrayListDble(rarraylist,itable,p_DData,ndata)

!<description>
    ! This subroutine copies the content of the list of the given
    ! table to the given double array
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)         :: rarraylist

    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)    :: itable
!</input>

!<inputoutput>
    ! double array
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p_DData
!</inputoutput>

!<output>
    ! OPTIONAL: number of data entries
    INTEGER, INTENT(OUT), OPTIONAL :: ndata
!</output>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,icount

    ! Check if table is valid
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) THEN
      CALL output_line('Invalid table number',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListDble')
      CALL sys_halt()
    END IF
    
    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListDble')
      CALL sys_halt()
    END IF

    icount = 0
    ipos = rarraylist%Ktable(ARRLST_HEAD,itable)
    DO
      icount = icount+1
      p_DData(icount) = rarraylist%DData(ipos)
      IF (ipos .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) EXIT
      ipos = rarraylist%Knext(ipos)
    END DO

    IF (PRESENT(ndata)) ndata=icount
  END SUBROUTINE arrlst_copyFromArrayListDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyFromArrayListDble_table(rarraylist,p_DData,p_Table)

!<description>
    ! This subroutine copies the content of the table and all lists
    ! to the given arrays.
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)                       :: rarraylist
!</input>

!<inputoutput>
    ! double array
    REAL(DP), DIMENSION(:), INTENT(INOUT)               :: p_DData

    ! table array
    INTEGER(PREC_TABLEIDX), DIMENSION(:), INTENT(INOUT) :: p_Table
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TABLEIDX) :: icount,itable,ntable
    INTEGER(PREC_ARRAYLISTIDX) :: ipos

    ! Check if table array is valid
    ntable=SIZE(p_Table)-1
    IF (ntable+1 < rarraylist%NTABLE) THEN
      CALL output_line('Invalid dimension of table array!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListDble_table')
      CALL sys_halt()
    END IF
    
    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListDble_table')
      CALL sys_halt()
    END IF

    icount=1
    DO itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%Ktable(ARRLST_HEAD,itable)
      DO WHILE (ipos .NE. ARRLST_NULL)
        p_DData(icount) = rarraylist%DData(ipos)
        icount = icount+1
        IF (ipos .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) EXIT
        ipos = rarraylist%Knext(ipos)
      END DO
    END DO
    p_Table(ntable+1)=icount+1
  END SUBROUTINE arrlst_copyFromArrayListDble_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyFromArrayListSngl(rarraylist,itable,p_SData,ndata)

!<description>
    ! This subroutine copies the content of the list of the given
    ! table to the given single array
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)         :: rarraylist
    
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)    :: itable
!</input>

!<inputoutput>
    ! double array
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: p_SData
!</inputoutput>

!<output>
    ! OPTIONAL: number of data entries
    INTEGER, INTENT(OUT), OPTIONAL :: ndata
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,icount

    ! Check if table is valid
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListSngl')
      CALL sys_halt()
    END IF

    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListSngl')
      CALL sys_halt()
    END IF

    icount = 0
    ipos = rarraylist%Ktable(ARRLST_HEAD,itable)
    DO
      icount = icount+1
      p_SData(icount) = rarraylist%SData(ipos)
      IF (ipos .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) EXIT
      ipos = rarraylist%Knext(ipos)
    END DO

    IF (PRESENT(ndata)) ndata=icount
  END SUBROUTINE arrlst_copyFromArrayListSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyFromArrayListSngl_table(rarraylist,p_SData,p_Table)

!<description>
    ! This subroutine copies the content of the table and all lists
    ! to the given arrays.   
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)                       :: rarraylist
!</input>

!<inputoutput>
    ! single array
    REAL(SP), DIMENSION(:), INTENT(INOUT)               :: p_SData

    ! table array
    INTEGER(PREC_TABLEIDX), DIMENSION(:), INTENT(INOUT) :: p_Table
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TABLEIDX) :: icount,itable,ntable
    INTEGER(PREC_ARRAYLISTIDX) :: ipos

    ! Check if table array is valid
    ntable=SIZE(p_Table)-1
    IF (ntable+1 < rarraylist%NTABLE) THEN
      CALL output_line('Invalid dimension of table array!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListSngl_table')
      CALL sys_halt()
    END IF
    
    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListSngl_table')
      CALL sys_halt()
    END IF

    icount=1
    DO itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%Ktable(ARRLST_HEAD,itable)
      DO WHILE(ipos .NE. ARRLST_NULL)
        p_SData(icount) = rarraylist%SData(ipos)
        icount = icount+1
        IF (ipos .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) EXIT
        ipos = rarraylist%Knext(ipos)
      END DO
    END DO
    p_Table(ntable+1)=icount+1
  END SUBROUTINE arrlst_copyFromArrayListSngl_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyFromArrayListInt(rarraylist,itable,p_IData,ndata)

!<description>
    ! This subroutine copies the content of the list of the given
    ! table to the given integer array.
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)        :: rarraylist

    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)   :: itable
!</input>

!<inputoutput>
    ! double array
    INTEGER, DIMENSION(:), INTENT(INOUT) :: p_IData
!</inputoutput>

!<output>
    ! OPTIONAL: number of data entries
    INTEGER, INTENT(OUT), OPTIONAL :: ndata
!</output>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,icount

    ! Check if table is valid
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListInt')
      CALL sys_halt()
    END IF

    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListInt')
      CALL sys_halt()
    END IF

    icount = 0
    ipos = rarraylist%Ktable(ARRLST_HEAD,itable)
    DO
      icount = icount+1
      p_IData(icount) = rarraylist%IData(ipos)
      IF (ipos .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) EXIT
      ipos = rarraylist%Knext(ipos)
    END DO

    IF (PRESENT(ndata)) ndata=icount
  END SUBROUTINE arrlst_copyFromArrayListInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyFromArrayListInt_table(rarraylist,p_IData,p_Table)

!<description>
    ! This subroutine copies the content of the table and all lists
    ! to the given arrays.   
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)                           :: rarraylist
!</input>

!<inputoutput>
    ! integer array
    INTEGER(PREC_ARRAYLISTIDX), DIMENSION(:), INTENT(INOUT) :: p_IData

    ! table array
    INTEGER(PREC_TABLEIDX), DIMENSION(:), INTENT(INOUT)     :: p_Table
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TABLEIDX) :: icount,itable,ntable
    INTEGER(PREC_ARRAYLISTIDX) :: ipos

    ! Check if table array is valid
    ntable=SIZE(p_Table)-1
    IF (ntable+1 < rarraylist%NTABLE) THEN
      CALL output_line('Invalid dimension of table array!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListInt_table')
      CALL sys_halt()
    END IF
    
    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyFromArrayListInt_table')
      CALL sys_halt()
    END IF

    icount=1
    DO itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%Ktable(ARRLST_HEAD,itable)
      DO WHILE(ipos .NE. ARRLST_NULL)
        p_IData(icount) = rarraylist%IData(ipos)
        icount=icount+1
        IF (ipos .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) EXIT
        ipos = rarraylist%Knext(ipos)
      END DO
    END DO
    p_Table(ntable+1)=icount
  END SUBROUTINE arrlst_copyFromArrayListInt_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyToArrayList(h_DataSrc,itable,rarraylist)

!<description>
    ! This subroutine copies the content of the given handle to the
    ! list associated with the given table
!</description>

!<input>
    ! handle to the data
    INTEGER, INTENT(IN)                :: h_DataSrc
     
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_DData
    REAL(SP), DIMENSION(:), POINTER :: p_SData
    INTEGER,  DIMENSION(:), POINTER :: p_IData
    
    ! Transform the content of h_Data to the list
    SELECT CASE (rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_DataSrc,p_DData)
      CALL arrlst_copyArrayList(p_DData,itable,rarraylist)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_DataSrc,p_SData)
      CALL arrlst_copyArrayList(p_SData,itable,rarraylist)

    CASE (ST_INT)
      CALL storage_getbase_int(h_DataSrc,p_IData)
      CALL arrlst_copyArrayList(p_IData,itable,rarraylist)
      stop

    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayList')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE arrlst_copyToArrayList

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyToArrayList_table(h_DataSrc,rarraylist,h_Table)

!<description>
    ! This subroutine copies the content of the given handle to the
    ! lists of the arraylist making use of the second handle to
    ! generate the table
!</description>

!<input>
    ! handle to the data
    INTEGER, INTENT(IN)              :: h_DataSrc

    ! handle to the table
    INTEGER, INTENT(IN)              :: h_Table
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TABLEIDX), DIMENSION(:), POINTER :: p_Table
    REAL(DP), DIMENSION(:), POINTER :: p_DData
    REAL(SP), DIMENSION(:), POINTER :: p_SData
    INTEGER,  DIMENSION(:), POINTER :: p_IData
    
    ! Set pointer to table
    CALL storage_getbase_int(h_Table,p_Table)

    ! Transform the content of h_Data to the list
    SELECT CASE (rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_DataSrc,p_DData)
      CALL arrlst_copyArrayListTable(p_DData,rarraylist,p_Table)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_DataSrc,p_SData)
      CALL arrlst_copyArrayListTable(p_SData,rarraylist,p_Table)

    CASE (ST_INT)
      CALL storage_getbase_int(h_DataSrc,p_IData)
      CALL arrlst_copyArrayListTable(p_IData,rarraylist,p_Table)

    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayList_table')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE arrlst_copyToArrayList_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyToArrayListDble(p_DDataSrc,itable,rarraylist)

!<description>
    ! This subroutine copies the content of the given double array to
    ! the list associated with the given table
!</description>

!<input>
    ! pointer to the data
    REAL(DP), DIMENSION(:), INTENT(IN) :: p_DDataSrc

    ! number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,kpos

    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListDble')
      CALL sys_halt()
    END IF
    
    DO ipos=1,SIZE(p_DDataSrc)
      CALL arrlst_appendToArrayList(rarraylist,itable,p_DDataSrc(ipos),kpos)
    END DO
  END SUBROUTINE arrlst_copyToArrayListDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyToArrayListDble_table(p_DDataSrc,rarraylist,p_Table)

!<description>
    ! This subroutine copies the content of the given double array to
    ! the lists of the arraylist making use of the second table array
!</description>

!<input>
    ! pointer to the data
    REAL(DP), DIMENSION(:), INTENT(IN)               :: p_DDataSrc

    ! pointer to the table
    INTEGER(PREC_TABLEIDX), DIMENSION(:), INTENT(IN) :: p_Table
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)                 :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_TABLEIDX)     :: itable,ntable
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,kpos

    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListDble_table')
      CALL sys_halt()
    END IF
    
    ntable=SIZE(p_Table)-1
    DO itable=1,ntable
      DO ipos=p_Table(itable),p_Table(itable+1)-1
        CALL arrlst_appendToArrayList(rarraylist,itable,p_DDataSrc(ipos),kpos)
      END DO
    END DO
  END SUBROUTINE arrlst_copyToArrayListDble_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyToArrayListSngl(p_SDataSrc,itable,rarraylist)

!<description>
    ! This subroutine copies the content of the given single array to
    ! the list associated with the given table
!</description>

!<input>
    ! pointer to the data
    REAL(SP), DIMENSION(:), INTENT(IN) :: p_SDataSrc

    ! number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,kpos

    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListSngl')
      CALL sys_halt()
    END IF
    
    DO ipos=1,SIZE(p_SDataSrc)
      CALL arrlst_appendToArrayList(rarraylist,itable,p_SDataSrc(ipos),kpos)
    END DO
  END SUBROUTINE arrlst_copyToArrayListSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyToArrayListSngl_table(p_SDataSrc,rarraylist,p_Table)

!<description>
    ! This subroutine copies the content of the given single array to
    ! the lists of the arraylist making use of the second table array
!</description>

!<input>
    ! pointer to the data
    REAL(SP), DIMENSION(:), INTENT(IN)               :: p_SDataSrc

    ! pointer to the table
    INTEGER(PREC_TABLEIDX), DIMENSION(:), INTENT(IN) :: p_Table
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)                 :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_TABLEIDX) :: itable,ntable
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,kpos

    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListSngl_table')
      CALL sys_halt()
    END IF

    ntable=SIZE(p_Table)-1
    DO itable=1,ntable
      DO ipos=p_Table(itable),p_Table(itable+1)-1
        CALL arrlst_appendToArrayList(rarraylist,itable,p_SDataSrc(ipos),kpos)
      END DO
    END DO
  END SUBROUTINE arrlst_copyToArrayListSngl_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyToArrayListInt(p_IDataSrc,itable,rarraylist)

!<description>
    ! This subroutine copies the content of the given integer array
    ! to the list associated with the given table
!</description>

!<input>
    ! pointer to the data
    INTEGER, DIMENSION(:), INTENT(IN)  :: p_IDataSrc

    ! number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,kpos

    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListInt')
      CALL sys_halt()
    END IF
    
    DO ipos=1,SIZE(p_IDataSrc)
      CALL arrlst_appendToArrayList(rarraylist,itable,p_IDataSrc(ipos),kpos)
    END DO
  END SUBROUTINE arrlst_copyToArrayListInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_copyToArrayListInt_table(p_IDataSrc,rarraylist,p_Table)

!<description>
    ! This subroutine copies the content of the given integer array
    ! to the lists of the arraylist making use of the second table array
!</description>

!<input>
    ! pointer to the data
    INTEGER, DIMENSION(:), INTENT(IN)                :: p_IDataSrc

    ! pointer to the table
    INTEGER(PREC_TABLEIDX), DIMENSION(:), INTENT(IN) :: p_Table
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)                 :: rarraylist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_TABLEIDX) :: itable,ntable
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,kpos

    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_copyToArrayListInt_table')
      CALL sys_halt()
    END IF
    
    ntable=SIZE(p_Table)-1
    DO itable=1,ntable
      DO ipos=p_Table(itable),p_Table(itable+1)-1
        CALL arrlst_appendToArrayList(rarraylist,itable,p_IDataSrc(ipos),kpos)
      END DO
    END DO
  END SUBROUTINE arrlst_copyToArrayListInt_table
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE arrlst_swapArrayList(rarraylist,itable,jtable)

!<description>
    ! This subroutine swaps two tables and the associated list in the
    ! given arraylist
!</description>
    
!<input>
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable,jtable
!</input>

!<inputoutput>
    ! arraylist
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ihead,itail,iitem,ina
    
    ! Swap
    ihead = rarraylist%Ktable(ARRLST_HEAD,itable)
    itail = rarraylist%Ktable(ARRLST_TAIL,itable)
    iitem = rarraylist%Ktable(ARRLST_ITEM,itable)
    ina   = rarraylist%Ktable(ARRLST_NA,  itable)
    
    rarraylist%Ktable(ARRLST_HEAD,itable) = rarraylist%Ktable(ARRLST_HEAD,jtable)
    rarraylist%Ktable(ARRLST_TAIL,itable) = rarraylist%Ktable(ARRLST_TAIL,jtable)
    rarraylist%Ktable(ARRLST_ITEM,itable) = rarraylist%Ktable(ARRLST_ITEM,jtable)
    rarraylist%Ktable(ARRLST_NA,  itable) = rarraylist%Ktable(ARRLST_NA,  jtable)

    rarraylist%Ktable(ARRLST_HEAD,jtable) = ihead
    rarraylist%Ktable(ARRLST_TAIL,jtable) = itail
    rarraylist%Ktable(ARRLST_ITEM,jtable) = iitem
    rarraylist%Ktable(ARRLST_NA,  jtable) = ina
  END SUBROUTINE arrlst_swapArrayList

  ! ***************************************************************************
  
!<function>
  
  PURE FUNCTION arrlst_getFirstInArrayList(rarraylist,itable) RESULT(ipos)

!<description>
    ! This function returns the position of the first list item in
    ! the given table. If the specified table does not exist than
    ! IPOS=-1 is returned.
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)      :: rarraylist

    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<result>
    ! position of first item
    INTEGER(PREC_ARRAYLISTIDX)         :: ipos
!</result>
!</function>

    IF (itable < 1 .OR. itable > rarraylist%NTABLE) THEN
      ipos = -1
      RETURN
    END IF
    
    ipos=rarraylist%Knext(rarraylist%Ktable(ARRLST_HEAD,itable))
  END FUNCTION arrlst_getFirstInArrayList

  ! ***************************************************************************
  
!<function>
  
  PURE FUNCTION arrlst_getLastInArrayList(rarraylist,itable) RESULT(ipos)

!<description>
    ! This function returns the position of the last list item in the
    ! given table. Of the specified table does not exist than IPOS=-1
    ! is returned.
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)      :: rarraylist

    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable
!</input>

!<result>
    ! position of last item
    INTEGER(PREC_ARRAYLISTIDX)         :: ipos
!</result>
!</function>

    IF (itable < 1 .OR. itable > rarraylist%NTABLE) THEN
      ipos = -1
      RETURN
    END IF
    
    ipos=rarraylist%Knext(rarraylist%Ktable(ARRLST_TAIL,itable))
  END FUNCTION arrlst_getLastInArrayList

  ! ***************************************************************************

!<function>

  FUNCTION arrlst_getNextInArrayList(rarraylist,itable,breset) RESULT(ipos)

!<description>
    ! This function returns the position of the next free list item
    ! for a given table and resets the list if required. If the
    ! specified table does not exist then IPOS=-1 is returned.
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable

    ! Reset list?
    LOGICAL, INTENT(IN)                :: breset
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>

!<result>
    ! position of last item
    INTEGER(PREC_ARRAYLISTIDX)         :: ipos
!</result>
!</function>
    
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) THEN
      ipos = -1
      RETURN
    END IF

    ! Should we reset the item pointer?
    IF (breset) THEN
      ipos = rarraylist%Ktable(ARRLST_HEAD,itable)
      rarraylist%Ktable(ARRLST_ITEM,itable) = ipos
      RETURN
    END IF

    ! Get next item and increase item pointer
    ipos = rarraylist%Knext(rarraylist%Ktable(ARRLST_ITEM,itable))
    rarraylist%Ktable(ARRLST_ITEM,itable) = ipos
  END FUNCTION arrlst_getNextInArrayList

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_prependToArrayListDble(rarraylist,itable,da,ipos)

!<description>
    ! This subroutine prepends a Double data to the list of the given
    ! table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! Data
    REAL(DP), INTENT(IN)                    :: da
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListDble')
      CALL sys_halt()
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListDble')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)
    
    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
   
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%DData(ipos)                = da
    ELSE
      rarraylist%Knext(ipos)                = rarraylist%Ktable(ARRLST_HEAD,itable)
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%DData(ipos)                = da
    END IF
  END SUBROUTINE arrlst_prependToArrayListDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_prependToArrayListSngl(rarraylist,itable,sa,ipos)

!<description>
    ! This subroutine prepends a Single data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! Data
    REAL(SP), INTENT(IN)                    :: sa
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListSngl')
      CALL sys_halt()
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListSngl')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
     rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%SData(ipos)                = sa
    ELSE
      rarraylist%Knext(ipos)                = rarraylist%Ktable(ARRLST_HEAD,itable)
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%SData(ipos)                = sa
    END IF
  END SUBROUTINE arrlst_prependToArrayListSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_prependToArrayListInt(rarraylist,itable,ia,ipos)

!<description>
    ! This subroutine prepends an Integer data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! Data
    INTEGER, INTENT(IN)                     :: ia
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListInt')
      CALL sys_halt()
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_prependToArrayListInt')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
     rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%IData(ipos)                = ia
    ELSE
      rarraylist%Knext(ipos)                = rarraylist%Ktable(ARRLST_HEAD,itable)
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%IData(ipos)                = ia
    END IF
  END SUBROUTINE arrlst_prependToArrayListInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_appendToArrayListDble(rarraylist,itable,da,ipos)

!<description>
    ! This subroutine appends a Double data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! Data
    REAL(DP), INTENT(IN)                    :: da
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the appended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListDble')
      CALL sys_halt()
    END IF

    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListDble')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)
    
    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%DData(ipos)                = da
    ELSE
      rarraylist%Knext(rarraylist%Ktable(ARRLST_TAIL,itable)) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable)                   = ipos
      rarraylist%Knext(ipos)                                  = ARRLST_NULL
      rarraylist%DData(ipos)                                  = da
    END IF
  END SUBROUTINE arrlst_appendToArrayListDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_appendToArrayListSngl(rarraylist,itable,sa,ipos)

!<description>
    ! This subroutine appends a Single data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable
    
    ! Data
    REAL(SP), INTENT(IN)                    :: sa
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the appended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListSngl')
      CALL sys_halt()
    END IF

    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListSngl')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%SData(ipos)                = sa
    ELSE
      rarraylist%Knext(rarraylist%Ktable(ARRLST_TAIL,itable)) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable)                   = ipos
      rarraylist%Knext(ipos)                                  = ARRLST_NULL
      rarraylist%SData(ipos)                                  = sa
    END IF
  END SUBROUTINE arrlst_appendToArrayListSngl
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_appendToArrayListInt(rarraylist,itable,ia,ipos)

!<description>
    ! This subroutine appends an Integer data to the list of a given
    ! table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! Data
    INTEGER, INTENT(IN)                     :: ia
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the appended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListInt')
      CALL sys_halt()
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_appendToArrayListInt')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%IData(ipos)                = ia
    ELSE
      rarraylist%Knext(rarraylist%Ktable(ARRLST_TAIL,itable)) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable)                   = ipos
      rarraylist%Knext(ipos)                                  = ARRLST_NULL
      rarraylist%IData(ipos)                                  = ia
    END IF
  END SUBROUTINE arrlst_appendToArrayListInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_insertIntoArrayListDble(rarraylist,itable,da,ipred,ipos)

!<description>
    ! This subroutine inserts a new Double data into the list of a
    ! given table AFTER the position ipred
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! Data
    REAL(DP), INTENT(IN)                    :: da

    ! Position of predecessor
    INTEGER(PREC_ARRAYLISTIDX), INTENT(IN)  :: ipred
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListDble')
      CALL sys_halt()
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListDble')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%DData(ipos)                = da
    ELSEIF (ipred .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) THEN
      rarraylist%Knext(ipred)               = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%DData(ipos)                = da
    ELSEIF (ipred < ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Knext(ipos)                = -ipred
      rarraylist%DData(ipos)                = da
    ELSE
      rarraylist%Knext(ipos)                = rarraylist%Knext(ipred)
      rarraylist%Knext(ipred)               = ipos
      rarraylist%DData(ipos)                = da
    END IF
  END SUBROUTINE arrlst_insertIntoArrayListDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_insertIntoArrayListSngl(rarraylist,itable,sa,ipred,ipos)

!<description>
    ! This subroutine inserts a new Single data into the list of a
    ! given table AFTER the position ipred
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! Data
    REAL(SP), INTENT(IN)                    :: sa

    ! Position of predecessor
    INTEGER(PREC_ARRAYLISTIDX), INTENT(IN)  :: ipred
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListSngl')
      CALL sys_halt()
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListSngl')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%SData(ipos)                = sa
    ELSEIF (ipred .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) THEN
      rarraylist%Knext(ipred)               = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%SData(ipos)                = sa
    ELSEIF (ipred < ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Knext(ipos)                = -ipred
      rarraylist%SData(ipos)                = sa
    ELSE
      rarraylist%Knext(ipos)                = rarraylist%Knext(ipred)
      rarraylist%Knext(ipred)               = ipos
      rarraylist%SData(ipos)                = sa
    END IF
  END SUBROUTINE arrlst_insertIntoArrayListSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_insertIntoArrayListInt(rarraylist,itable,ia,ipred,ipos)

!<description>
    ! This subroutine inserts a new Integer data into the list of a
    ! given table AFTER the position ipred
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! Data
    INTEGER, INTENT(IN)                     :: ia

    ! Position of predecessor
    INTEGER(PREC_ARRAYLISTIDX), INTENT(IN)  :: ipred
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)        :: rarraylist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListInt')
      CALL sys_halt()
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      CALL output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_insertIntoArrayListInt')
      CALL sys_halt()
    END IF
    IF (rarraylist%NTABLE < itable) CALL arrlst_createArrayList_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%Ktable(ARRLST_NA,itable) =  rarraylist%Ktable(ARRLST_NA,itable)+1
    ipos = rarraylist%Knext(ARRLST_FREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL arrlst_resizeArrayList(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(ARRLST_FREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(ARRLST_FREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(ARRLST_HEAD,itable) .EQ. ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%IData(ipos)                = ia
    ELSEIF (ipred .EQ. rarraylist%Ktable(ARRLST_TAIL,itable)) THEN
      rarraylist%Knext(ipred)               = ipos
      rarraylist%Ktable(ARRLST_TAIL,itable) = ipos
      rarraylist%Knext(ipos)                = ARRLST_NULL
      rarraylist%IData(ipos)                = ia
    ELSEIF (ipred < ARRLST_NULL) THEN
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Knext(ipos)                = -ipred
      rarraylist%IData(ipos)                = ia
    ELSE
      rarraylist%Knext(ipos)                = rarraylist%Knext(ipred)
      rarraylist%Knext(ipred)               = ipos
      rarraylist%IData(ipos)                = ia
    END IF
  END SUBROUTINE arrlst_insertIntoArrayListInt
  
  ! ***************************************************************************
  
!<function>

  FUNCTION arrlst_deleteFromArrayListDble(rarraylist,itable,da) RESULT(f)

!<description>
    ! This function deletes a Double data from the arraylist
    ! associated with the given table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable

    ! Data
    REAL(DP), INTENT(IN)               :: da
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / ARRAYLIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipred,ipos

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_deleteFromArrayListDble')
      CALL sys_halt()
    END IF

    ! Search for data
    f=arrlst_searchInArrayList(rarraylist,itable,da,ipred)
    IF (f .EQ. ARRAYLIST_NOT_FOUND) RETURN
    
    ! Delete data
    rarraylist%NA = rarraylist%NA-1
    rarraylist%Ktable(ARRLST_NA,itable) = rarraylist%Ktable(ARRLST_NA,itable)-1

    ! Are we first entry in list?
    IF (ipred < 0) THEN

      ! Get position
      ipred = -ipred
      ipos  = rarraylist%Knext(ipred)
      
      ! Update free position
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Knext(ipred)               = rarraylist%Knext(ARRLST_FREE)
      rarraylist%Knext(ARRLST_FREE)         = -ipred

    ELSE

      ! Get position
      ipos = rarraylist%Knext(ipred)
      IF (rarraylist%Knext(ipred) .EQ. rarraylist%Ktable(ARRLST_TAIL,itable))&
          rarraylist%Ktable(ARRLST_TAIL,itable) = ipred
      
      ! Update free position
      rarraylist%Knext(ipred)       = rarraylist%Knext(ipos)
      rarraylist%Knext(ipos)        = rarraylist%Knext(ARRLST_FREE)
      rarraylist%Knext(ARRLST_FREE) = -ipos
    END IF
  END FUNCTION arrlst_deleteFromArrayListDble
  
  ! ***************************************************************************
  
!<function>

  FUNCTION arrlst_deleteFromArrayListSngl(rarraylist,itable,sa) RESULT(f)

!<description>
    ! This function deletes a Single data from the arraylist
    ! associated with the given table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable

    ! Data
    REAL(SP), INTENT(IN)               :: sa
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / ARRAYLIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipred,ipos

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_deleteFromArrayListSngl')
      CALL sys_halt()
    END IF

    ! Search for data
    f=arrlst_searchInArrayList(rarraylist,itable,sa,ipred)
    IF (f .EQ. ARRAYLIST_NOT_FOUND) RETURN

    ! Delete data
    rarraylist%NA = rarraylist%NA-1
    rarraylist%Ktable(ARRLST_NA,itable) = rarraylist%Ktable(ARRLST_NA,itable)-1
    
    ! Are we first entry in list?
    IF (ipred < 0) THEN
      
      ! Get position
      ipred = -ipred
      ipos  = rarraylist%Knext(ipred)
      
      ! Update free position
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Knext(ipred)         = rarraylist%Knext(ARRLST_FREE)
      rarraylist%Knext(ARRLST_FREE)         = -ipred
      
    ELSE
      
      ! Get position
      ipos = rarraylist%Knext(ipred)
      IF (rarraylist%Knext(ipred) .EQ. rarraylist%Ktable(ARRLST_TAIL,itable))&
          rarraylist%Ktable(ARRLST_TAIL,itable)=ipred
      
      ! Update free position
      rarraylist%Knext(ipred) = rarraylist%Knext(ipos)
      rarraylist%Knext(ipos)  = rarraylist%Knext(ARRLST_FREE)
      rarraylist%Knext(ARRLST_FREE) = -ipos
    END IF
  END FUNCTION arrlst_deleteFromArrayListSngl

  ! ***************************************************************************
  
!<function>

  FUNCTION arrlst_deleteFromArrayListInt(rarraylist,itable,ia) RESULT(f)

!<description>
    ! This function deletes an Integer data from the arraylist
    ! associated with the given table.
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN) :: itable

    ! Data
    INTEGER, INTENT(IN)                :: ia
!</input>

!<inputoutput>
    ! list
    TYPE(t_arraylist), INTENT(INOUT)   :: rarraylist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / ARRAYLIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipred,ipos

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_deleteFromArrayListInt')
      CALL sys_halt()
    END IF

    ! Search for data
    f=arrlst_searchInArrayList(rarraylist,itable,ia,ipred)
    IF (f .EQ. ARRAYLIST_NOT_FOUND) RETURN

    ! Delete data
    rarraylist%NA = rarraylist%NA-1
    rarraylist%Ktable(ARRLST_NA,itable) = rarraylist%Ktable(ARRLST_NA,itable)-1

    ! Are we first entry in list?
    IF (ipred < 0) THEN

      ! Get position
      ipred = -ipred
      ipos  = rarraylist%Knext(ipred)

      ! Update free position
      rarraylist%Ktable(ARRLST_HEAD,itable) = ipos
      rarraylist%Knext(ipred)               = rarraylist%Knext(ARRLST_FREE)
      rarraylist%Knext(ARRLST_FREE)         = -ipred
      
    ELSE

      ! Get position
      ipos = rarraylist%Knext(ipred)

      ! Check if last entry should be deleted
      IF (rarraylist%Knext(ipred) .EQ. rarraylist%Ktable(ARRLST_TAIL,itable))&
          rarraylist%Ktable(ARRLST_TAIL,itable)=ipred
      
      ! Update free position
      rarraylist%Knext(ipred) = rarraylist%Knext(ipos)
      rarraylist%Knext(ipos)  = rarraylist%Knext(ARRLST_FREE)
      rarraylist%Knext(ARRLST_FREE) = -ipos
    END IF
  END FUNCTION arrlst_deleteFromArrayListInt

  ! ***************************************************************************
  
!<function>

  FUNCTION arrlst_searchInArrayListDble(rarraylist,itable,da,ipred) RESULT(f)

!<description>
    ! This function searches for a given Double data in the list of a
    ! given table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! list
    TYPE(t_arraylist), INTENT(IN)           :: rarraylist

    ! Data
    REAL(DP), INTENT(IN)                    :: da
!</input>

!<output>
    ! Position of the predecessor of the found item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOUT_FOUND / ARRAYLIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ihead,itail,inext

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_searchInArrayListDble')
      CALL sys_halt()
    END IF

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) RETURN
    ipred=-rarraylist%Ktable(ARRLST_HEAD,itable)
    
    ! Check if list is empty
    IF (ipred .EQ. ARRLST_NULL) RETURN

    ! Initialization
    ihead=rarraylist%Ktable(ARRLST_HEAD,itable)
    itail=rarraylist%Ktable(ARRLST_TAIL,itable)

    ! What kind of ordering are we
    SELECT CASE(rarraylist%cordering)
    CASE (ARRAYLIST_UNORDERED)

      ! Check first item separately
      IF (rarraylist%DData(ihead) .EQ. da) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%DData(inext) .EQ. da) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (inext .EQ. itail) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_INCREASING)

      ! Check first item separately
      IF (rarraylist%DData(ihead) .EQ. da) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%DData(ihead) > da) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%DData(inext) .EQ. da) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%DData(inext) > da) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_DECREASING)

      ! Check first item separately
      IF (rarraylist%DData(ihead) .EQ. da) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%DData(ihead) < da) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)     
        IF (rarraylist%DData(inext) .EQ. da) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%DData(inext) < da) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_CSR7)

      ! Check first item separately
      IF (rarraylist%DData(ihead) .EQ. da) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(ABS(itable-da) .LE. SYS_EPSREAL) THEN
        RETURN
      END IF
      ipred=-ipred
      
      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%DData(inext) .EQ. da) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%DData(inext) > da) EXIT
        IF (rarraylist%Knext(ipred) .EQ. itail) THEN
          ipred=rarraylist%Knext(ipred); EXIT
        END IF
        ipred=rarraylist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION arrlst_searchInArrayListDble
  
  ! ***************************************************************************
  
!<function>

  FUNCTION arrlst_searchInArrayListSngl(rarraylist,itable,sa,ipred) RESULT(f)

!<description>
    ! This function searches for a given Single data in the list of a
    ! given table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! list
    TYPE(t_arraylist), INTENT(IN)           :: rarraylist

    ! Data
    REAL(SP), INTENT(IN)                    :: sa
!</input>

!<output>
    ! Position of the predecessor of the found item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipred
!</output>

!<result>
    ! Result of the search ARRAYLIST_NOT_FOUND / ARRAYLIST_FOUND
    INTEGER :: f
!</result>
!</function>
    
    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ihead,itail,inext

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_searchInArrayListSngl')
      CALL sys_halt()
    END IF

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) RETURN
    ipred=-rarraylist%Ktable(ARRLST_HEAD,itable)

    ! Check if list is empty
    IF (ipred .EQ. ARRLST_NULL) RETURN

    ! Initialization
    ihead=rarraylist%Ktable(ARRLST_HEAD,itable)
    itail=rarraylist%Ktable(ARRLST_TAIL,itable)

    ! What kind of ordering are we
    SELECT CASE(rarraylist%cordering)
    CASE (ARRAYLIST_UNORDERED)

      ! Check first item separately
      IF (rarraylist%SData(ihead) .EQ. sa) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%SData(inext) .EQ. sa) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (inext .EQ. itail) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_INCREASING)

      ! Check first item separately
      IF (rarraylist%SData(ihead) .EQ. sa) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%SData(ihead) > sa) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%SData(inext) .EQ. sa) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%SData(inext) > sa) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_DECREASING)

      ! Check first item separately
      IF (rarraylist%SData(ihead) .EQ. sa) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%SData(ihead) < sa) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%SData(inext) .EQ. sa) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%SData(inext) < sa) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_CSR7)

      ! Check first item separately
      IF (rarraylist%SData(ihead) .EQ. sa) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(ABS(itable-sa) .LE. SYS_EPSREAL) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%SData(inext) .EQ. sa) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%SData(inext) > sa) EXIT
        IF (rarraylist%Knext(ipred) .EQ. itail) THEN
          ipred=rarraylist%Knext(ipred); EXIT
        END IF
        ipred=rarraylist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION arrlst_searchInArrayListSngl

  ! ***************************************************************************
  
!<function>

  FUNCTION arrlst_searchInArrayListInt(rarraylist,itable,ia,ipred) RESULT(f)

!<description>
    ! This function searches for a given Integer data in the list of
    ! a given table
!</description>

!<input>
    ! Number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN)      :: itable

    ! list
    TYPE(t_arraylist), INTENT(IN)           :: rarraylist

    ! Data
    INTEGER, INTENT(IN)                     :: ia
!</input>

!<output>
    ! Position of the predecessor of the found item
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT) :: ipred
!</output>

!<result>
    ! Result of the search ARRAYLIST_NOT_FOUND / ARRAYLIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ihead,itail,inext

    ! Check if list format is ok
    IF (rarraylist%carraylistFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_searchInArrayListInt')
      CALL sys_halt()
    END IF

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) RETURN
    ipred=-rarraylist%Ktable(ARRLST_HEAD,itable)
    
    ! Check if list is empty
    IF (ipred .EQ. ARRLST_NULL) RETURN

    ! Initialization
    ihead=rarraylist%Ktable(ARRLST_HEAD,itable)
    itail=rarraylist%Ktable(ARRLST_TAIL,itable)
        
    ! What kind of ordering are we
    SELECT CASE(rarraylist%cordering)
    CASE (ARRAYLIST_UNORDERED)

      ! Check first item separately
      IF (rarraylist%IData(ihead) .EQ. ia) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%IData(inext) .EQ. ia) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF

        IF (inext .EQ. itail) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_INCREASING)

      ! Check first item separately
      IF (rarraylist%IData(ihead) .EQ. ia) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%IData(ihead) > ia) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%IData(inext) .EQ. ia) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%IData(inext) > ia) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_DECREASING)

      ! Check first item separately
      IF (rarraylist%IData(ihead) .EQ. ia) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%IData(ihead) < ia) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%IData(inext) .EQ. ia) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%IData(inext) < ia) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_CSR7)

      ! Check first item separately
      IF (rarraylist%IData(ihead) .EQ. ia) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(itable .EQ. ia) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%IData(inext) .EQ. ia) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%IData(inext) > ia) EXIT
        IF (rarraylist%Knext(ipred) .EQ. itail) THEN
          ipred=rarraylist%Knext(ipred); EXIT
        END IF
        ipred=rarraylist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION arrlst_searchInArrayListInt

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE arrlst_printArrayList(rarraylist,itable)

!<description>
    ! This subroutine prints the content of the list
!</description>

!<input>
    ! list
    TYPE(t_arraylist), INTENT(IN)                :: rarraylist

    ! OPTIONAL: number of table
    INTEGER(PREC_TABLEIDX), INTENT(IN), OPTIONAL :: itable
!</input>
!</subroutine>

    ! local variable
    INTEGER(PREC_TABLEIDX)     :: iitable,itable1,itable2
    INTEGER(PREC_ARRAYLISTIDX) :: ipos,itail
    
    IF (PRESENT(itable)) THEN
      itable1=itable; itable2=itable
    ELSE
      itable1=1; itable2=rarraylist%NTABLE
    END IF

    SELECT CASE (rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)

      CALL output_line('Array list:')
      CALL output_line('-----------')

      ! Loop over all selected tables
      DO iitable=itable1,itable2
        CALL output_line('Table number'//TRIM(sys_siL(iitable,15)))
        CALL output_line('----------------------------------------')

        ipos  = rarraylist%Ktable(ARRLST_HEAD,iitable)
        IF (ipos .EQ. ARRLST_NULL) CYCLE
        itail = rarraylist%Ktable(ARRLST_TAIL,iitable)
        
        DO
          WRITE(*,FMT='(A,",")',ADVANCE='NO') sys_sdL(rarraylist%DData(ipos),8)
          IF (ipos .EQ. itail) EXIT
          ipos = rarraylist%Knext(ipos)
        END DO
      END DO
      
    CASE (ST_SINGLE)

      CALL output_line('Array list:')
      CALL output_line('-----------')

      ! Loop over all selected tables
      DO iitable=itable1,itable2
        CALL output_line('Table number'//TRIM(sys_siL(iitable,15)))
        CALL output_line('----------------------------------------')
        
        ipos = rarraylist%Ktable(ARRLST_HEAD,iitable)
        IF (ipos .EQ. ARRLST_NULL) CYCLE
        itail = rarraylist%Ktable(ARRLST_TAIL,iitable)
        
        DO
          WRITE(*,FMT='(A,",")',ADVANCE='NO') sys_sdL(REAL(rarraylist%SData(ipos),DP),8)
          IF (ipos .EQ. itail) EXIT
          ipos = rarraylist%Knext(ipos)
        END DO
      END DO
      
    CASE (ST_INT)
      
      CALL output_line('Array list:')
      CALL output_line('-----------')
      
      ! Loop over all selected tables
      DO iitable=itable1,itable2
        CALL output_line('Table number'//TRIM(sys_siL(iitable,15)))
        CALL output_line('----------------------------------------')

        ipos = rarraylist%Ktable(ARRLST_HEAD,iitable)
        IF (ipos .EQ. ARRLST_NULL) CYCLE
        itail = rarraylist%Ktable(ARRLST_TAIL,iitable)
        
        DO
          WRITE(*,FMT='(A,",")',ADVANCE='NO') sys_siL(rarraylist%IData(ipos),15)
          IF (ipos .EQ. itail) EXIT
          ipos = rarraylist%Knext(ipos)
        END DO
        WRITE(*,*)
      END DO
      
    CASE DEFAULT
      CALL output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'arrlst_printArrayList')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE arrlst_printArrayList

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE arrlst_infoArrayList(rarraylist)

!<description>
    ! This subroutine prints information about the arraylist
!</description>

!<input>
    ! arraylist
    TYPE(t_arraylist), INTENT(IN) :: rarraylist
!</input>
!</subroutine>

    CALL output_line('Arraylist:')
    CALL output_line('----------')
    CALL output_line('NA:       '//TRIM(sys_siL(rarraylist%NA,15)))
    CALL output_line('NNA:      '//TRIM(sys_siL(rarraylist%NNA,15)))
    CALL output_line('NNA0:     '//TRIM(sys_siL(rarraylist%NNA0,15)))
    CALL output_line('NTABLE:   '//TRIM(sys_siL(rarraylist%NTABLE,15)))
    CALL output_line('NNTABLE:  '//TRIM(sys_siL(rarraylist%NNTABLE,15)))
    CALL output_line('NNTABLE0: '//TRIM(sys_siL(rarraylist%NNTABLE0,15)))
    CALL output_line('NRESIZE:  '//TRIM(sys_siL(rarraylist%NRESIZE,15)))
    CALL output_line('dfactor:  '//TRIM(sys_sdL(rarraylist%dfactor,2)))
    CALL output_line('h_Ktable: '//TRIM(sys_siL(rarraylist%h_Ktable,15)))
    CALL output_line('h_Knext:  '//TRIM(sys_siL(rarraylist%h_Knext,15)))
    CALL output_line('h_Data:   '//TRIM(sys_siL(rarraylist%h_Data,15)))

    CALL output_line('Current data  memory usage: '//&
        TRIM(sys_sdL(100*rarraylist%NA/REAL(rarraylist%NNA,DP),2))//'%')
    CALL output_line('Current table memory usage: '//&
        TRIM(sys_sdL(100*rarraylist%NTABLE/REAL(rarraylist%NNTABLE,DP),2))//'%')

    CALL output_line('Total   data  memory usage: '//&
        TRIM(sys_sdL(100*rarraylist%NNA/REAL(rarraylist%NNA0,DP),2))//'%')
    CALL output_line('Total   table memory usage: '//&
        TRIM(sys_sdL(100*rarraylist%NNTABLE/REAL(rarraylist%NNTABLE0,DP),2))//'%')
    CALL output_line('')

  END SUBROUTINE arrlst_infoArrayList
END MODULE arraylist
