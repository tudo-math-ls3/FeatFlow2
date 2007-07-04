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
!# 2.) arrlst_releaseArrayList = t_arraylist_release /
!#                               t_arraylist_release_table
!#     -> Release an existing arraylist
!#
!# 3.) arrlst_resizeArrayList = t_arraylist_resize
!#     -> Reallocate memory for an existing arraylist
!#
!# 4.) arrlst_copyArrayList = t_arraylist_copyfrom /
!#                            t_arraylist_copyfromDble /
!#                            t_arraylist_copyfromSngl /
!#                            t_arraylist_copyfromInt /
!#                            t_arraylist_copyto /
!#                            t_arraylist_copytoDble /
!#                            t_arraylist_copytoSngl /
!#                            t_arraylist_copytoInt
!#     -> Copy data to/from an arraylist for a given table
!#
!# 5.) arrlst_copyArrayListTable = t_arraylist_copyfrom_table /
!#                                 t_arraylist_copyfromDble_table /
!#                                 t_arraylist_copyfromSngl_table /
!#                                 t_arraylist_copyfromInt_table /
!#                                 t_arraylist_copyto_table /
!#                                 t_arraylist_copytoDble_table /
!#                                 t_arraylist_copytoSngl_table /
!#                                 t_arraylist_copytoInt_table
!#     -> Copy data to/from an arraylist for a complete table
!#        structure
!#
!# 6.) arrlst_swapArrayList
!#     -> Swap two lists in the table
!#
!# 7.) arrlst_getFirstInArraylist
!#     -> Get the position of the first item in arraylist
!#
!# 8.) arrlst_getLastInArraylist
!#     -> Get the position of the last item in arraylist
!#
!# 9.) arrlst_getNextInArraylist
!#     -> Get the position of the next item in arraylist
!#
!# 10.) arrlst_prependToArraylist = t_arraylist_prependDble /
!#                                  t_arraylist prependSngl /
!#                                  t_arraylist_prependInt
!#     -> Prepend data to arraylist
!#
!# 11.) arrlst_appendToArrayist = t_arraylist_appendDble /
!#                                t_arraylist_appendSngl /
!#                                t_arraylist_appendInt
!#      -> Append data to arraylist
!#
!# 12.) arrlst_insertIntoArraylist = t_arraylist_insertDble /
!#                                   t_arraylist_insertSngl /
!#                                   t_arraylist_insertInt
!#      -> Insert data into arraylist
!#
!# 13.) arrlst_deleteFromArraylist = t_arraylist_deleteDble /
!#                                   t_arraylist_deleteSngl /
!#                                   t_arraylist_deleteInt
!#      -> Delete data from arraylist
!#
!# 14.) arrlst_searchInArraylist = t_arraylist_searchDble /
!#                                 t_arraylist_searchSngl /
!#                                 t_arraylist_searchInt
!#      -> Search for data in arraylist
!#
!# 15.) arrlst_printArraylist
!#      -> Print content of arraylist
!#
!# 16.) arrlst_infoArraylist
!#      -> Prit info about arraylist
!#
!# </purpose>
!##############################################################################

MODULE arraylist
  USE fsystem
  USE storage
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_arraylist
  PUBLIC :: arrlst_createArraylist
  PUBLIC :: arrlst_releaseArraylist
  PUBLIC :: arrlst_resizeArraylist
  PUBLIC :: arrlst_copyArraylist
  PUBLIC :: arrlst_copyArraylistTable
  PUBLIC :: arrlst_swapArraylist
  PUBLIC :: arrlst_getFirstInArraylist
  PUBLIC :: arrlst_getLastInArraylist
  PUBLIC :: arrlst_getNextInArraylist
  PUBLIC :: arrlst_prependToArraylist
  PUBLIC :: arrlst_appendToArraylist
  PUBLIC :: arrlst_insertIntoArraylist
  PUBLIC :: arrlst_deleteFromArraylist
  PUBLIC :: arrlst_searchInArraylist
  PUBLIC :: arrlst_printArraylist
  PUBLIC :: arrlst_infoArraylist

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
  INTEGER, PARAMETER, PUBLIC :: ANULL =  0

  ! Tag for next free position in storage
  INTEGER, PARAMETER :: LFREE = 0

  ! Tag for head of arraylist
  INTEGER, PARAMETER :: LHEAD = 1

  ! Tag for tail of arraylist
  INTEGER, PARAMETER :: LTAIL = 2

  ! Tag for last item stored in arraylist
  INTEGER, PARAMETER :: LITEM = 3

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
    
    ! Arraylist structure
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_ARRAYLISTIDX), DIMENSION(:), POINTER :: Knext => NULL()

    ! Arraylist data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(DP), DIMENSION(:), POINTER :: DData => NULL()

    ! Arraylist data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(SP), DIMENSION(:), POINTER :: SData => NULL()

    ! Arraylist data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_ARRAYLISTIDX), DIMENSION(:), POINTER :: IData => NULL()
  END TYPE t_arraylist
  
!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  INTERFACE arrlst_createArraylist
    MODULE PROCEDURE t_arraylist_create
  END INTERFACE

  INTERFACE create_table ! internal use
    MODULE PROCEDURE t_arraylist_create_table
  END INTERFACE
  
  INTERFACE arrlst_releaseArraylist
    MODULE PROCEDURE t_arraylist_release
    MODULE PROCEDURE t_arraylist_release_table
  END INTERFACE
  
  INTERFACE arrlst_resizeArraylist
    MODULE PROCEDURE t_arraylist_resize
  END INTERFACE
  INTERFACE resize   ! internal use
    MODULE PROCEDURE t_arraylist_resize
  END INTERFACE
  INTERFACE resize_table ! internal use
    MODULE PROCEDURE t_arraylist_resize_table
  END INTERFACE
  
  INTERFACE arrlst_copyArraylist
    MODULE PROCEDURE t_arraylist_copyfrom
    MODULE PROCEDURE t_arraylist_copyfromDble
    MODULE PROCEDURE t_arraylist_copyfromSngl
    MODULE PROCEDURE t_arraylist_copyfromInt
    MODULE PROCEDURE t_arraylist_copyto
    MODULE PROCEDURE t_arraylist_copytoDble
    MODULE PROCEDURE t_arraylist_copytoSngl
    MODULE PROCEDURE t_arraylist_copytoInt
  END INTERFACE

  INTERFACE arrlst_copyArraylistTable
    MODULE PROCEDURE t_arraylist_copyfrom_table
    MODULE PROCEDURE t_arraylist_copyfromDble_table
    MODULE PROCEDURE t_arraylist_copyfromSngl_table
    MODULE PROCEDURE t_arraylist_copyfromInt_table
    MODULE PROCEDURE t_arraylist_copyto_table
    MODULE PROCEDURE t_arraylist_copytoDble_table
    MODULE PROCEDURE t_arraylist_copytoSngl_table
    MODULE PROCEDURE t_arraylist_copytoInt_table
  END INTERFACE
  
  INTERFACE arrlst_swapArraylist
    MODULE PROCEDURE t_arraylist_swap
  END INTERFACE
  
  INTERFACE arrlst_getFirstInArraylist
    MODULE PROCEDURE t_arraylist_first
  END INTERFACE
  
  INTERFACE arrlst_getLastInArraylist
    MODULE PROCEDURE t_arraylist_last
  END INTERFACE
  
  INTERFACE arrlst_getNextInArraylist
    MODULE PROCEDURE t_arraylist_next
  END INTERFACE
  
  INTERFACE arrlst_prependToArraylist
    MODULE PROCEDURE t_arraylist_prependDble
    MODULE PROCEDURE t_arraylist_prependSngl
    MODULE PROCEDURE t_arraylist_prependInt
  END INTERFACE
  INTERFACE prepend   ! internal use
    MODULE PROCEDURE t_arraylist_prependDble
    MODULE PROCEDURE t_arraylist_prependSngl
    MODULE PROCEDURE t_arraylist_prependInt
  END INTERFACE
  
  INTERFACE arrlst_appendToArraylist
    MODULE PROCEDURE t_arraylist_appendDble
    MODULE PROCEDURE t_arraylist_appendSngl
    MODULE PROCEDURE t_arraylist_appendInt
  END INTERFACE
  INTERFACE append   ! internal use
    MODULE PROCEDURE t_arraylist_appendDble
    MODULE PROCEDURE t_arraylist_appendSngl
    MODULE PROCEDURE t_arraylist_appendInt
  END INTERFACE
  
  INTERFACE arrlst_insertIntoArraylist
    MODULE PROCEDURE t_arraylist_insertDble
    MODULE PROCEDURE t_arraylist_insertSngl
    MODULE PROCEDURE t_arraylist_insertInt
  END INTERFACE
  INTERFACE insert   ! internal use
    MODULE PROCEDURE t_arraylist_insertDble
    MODULE PROCEDURE t_arraylist_insertSngl
    MODULE PROCEDURE t_arraylist_insertInt
  END INTERFACE

  INTERFACE arrlst_deleteFromArraylist
    MODULE PROCEDURE t_arraylist_deleteDble
    MODULE PROCEDURE t_arraylist_deleteSngl
    MODULE PROCEDURE t_arraylist_deleteInt
  END INTERFACE
  INTERFACE delete   ! internal use
    MODULE PROCEDURE t_arraylist_deleteDble
    MODULE PROCEDURE t_arraylist_deleteSngl
    MODULE PROCEDURE t_arraylist_deleteInt
  END INTERFACE
  
  INTERFACE arrlst_searchInArraylist
    MODULE PROCEDURE t_arraylist_searchDble
    MODULE PROCEDURE t_arraylist_searchSngl
    MODULE PROCEDURE t_arraylist_searchInt
  END INTERFACE
  INTERFACE search   ! internal use
    MODULE PROCEDURE t_arraylist_searchDble
    MODULE PROCEDURE t_arraylist_searchSngl
    MODULE PROCEDURE t_arraylist_searchInt
  END INTERFACE
  
  INTERFACE arrlst_printArraylist
    MODULE PROCEDURE t_arraylist_print
  END INTERFACE

  INTERFACE arrlst_infoArraylist
    MODULE PROCEDURE t_arraylist_info
  END INTERFACE
  
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

CONTAINS
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE t_arraylist_create(rarraylist,nntable,nna,carraylistFormat,cordering,dfactor)

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
    Isize=(/3,nntable/)
    CALL storage_new('t_arraylist_create','Ktable',Isize,&
        ST_INT,rarraylist%h_Ktable,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int2D(rarraylist%h_Ktable,rarraylist%Ktable)
    
    CALL storage_new('t_arraylist_create','Knext',LFREE,nna,ST_INT,&
        rarraylist%h_Knext,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(rarraylist%h_Knext,rarraylist%Knext)

    SELECT CASE(rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_new('t_arraylist_create','Data',nna,ST_DOUBLE,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double(rarraylist%h_Data,rarraylist%DData)
      
    CASE (ST_SINGLE)
      CALL storage_new('t_arraylist_create','Data',nna,ST_SINGLE,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_single(rarraylist%h_Data,rarraylist%SData)
      
    CASE (ST_INT)
      CALL storage_new('t_arraylist_create','Data',nna,ST_INT,&
          rarraylist%h_Data,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int(rarraylist%h_Data,rarraylist%IData)
      
    CASE DEFAULT
      PRINT *, 't_arraylist_create: Unsupported data format!'
      STOP
    END SELECT
    
    ! Initialize list structures
    rarraylist%Knext(LFREE)    = 1
  END SUBROUTINE t_arraylist_create
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_create_table(rarraylist,itable)

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
      PRINT *, "t_arraylist_create_table: Invalid table number"
      STOP
    END IF

    ! Resize tables if required
    IF (itable > rarraylist%NNTABLE) CALL resize_table(rarraylist,&
        CEILING(itable*rarraylist%dfactor))
    
    ! Initialize structures
    rarraylist%Ktable(LHEAD:LITEM,rarraylist%NTABLE+1:itable) = ANULL

    ! Set new table size
    rarraylist%NTABLE=itable
  END SUBROUTINE t_arraylist_create_table

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_arraylist_release(rarraylist)

!<description>
    ! This subroutine releases an existing arraylist
!</description>

!<inputoutput>
    TYPE(t_arraylist), INTENT(INOUT) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Release memory
    IF (rarraylist%h_Ktable /= ST_NOHANDLE) CALL storage_free(rarraylist%h_Ktable)
    IF (rarraylist%h_Knext /= ST_NOHANDLE)  CALL storage_free(rarraylist%h_Knext)
    IF (rarraylist%h_Data /= ST_NOHANDLE)   CALL storage_free(rarraylist%h_Data)
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
  END SUBROUTINE t_arraylist_release

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_release_table(rarraylist,itable)

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
      PRINT *, "t_arraylist_release_table: Invalid table number!"
      STOP
    END IF

    ! Reset table
    rarraylist%Ktable(LHEAD:LITEM,itable) = ANULL
  END SUBROUTINE t_arraylist_release_table

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_arraylist_resize(rarraylist,nna)

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
    rarraylist%NNA=nna
    rarraylist%NRESIZE=rarraylist%NRESIZE+1

    CALL storage_realloc('t_arraylist_resize',LFREE,nna,&
        rarraylist%h_Knext,ST_NEWBLOCK_NOINIT,.TRUE.)
    CALL storage_realloc('t_arraylist_resize',nna,&
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
      PRINT *, 't_arraylist_resize: Unsupported data format!'
      STOP
    END SELECT
  END SUBROUTINE t_arraylist_resize

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_resize_table(rarraylist,nntable)

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
    rarraylist%NNTABLE=nntable
    rarraylist%NRESIZE=rarraylist%NRESIZE+1

    CALL storage_realloc('t_arraylist_resize_table',nntable,&
        rarraylist%h_Ktable,ST_NEWBLOCK_NOINIT,.TRUE.)
    CALL storage_getbase_int2D(rarraylist%h_Ktable,rarraylist%Ktable)

    rarraylist%Ktable(LHEAD,rarraylist%NTABLE+1:) = ANULL
    rarraylist%Ktable(LTAIL,rarraylist%NTABLE+1:) = ANULL
    rarraylist%Ktable(LITEM,rarraylist%NTABLE+1:) = LHEAD
  END SUBROUTINE t_arraylist_resize_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyfrom(rarraylist,itable,h_Data)

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
    IF (h_Data == ST_NOHANDLE) THEN
      CALL storage_new('t_arraylist_copyfrom','Data',rarraylist%na,&
          rarraylist%carraylistFormat,h_Data,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Data,isize)
      IF (isize < rarraylist%na) THEN
        CALL storage_realloc('t_arraylist_copyfrom',rarraylist%na,&
            h_Data,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF
    
    ! What kind of data are we?
    SELECT CASE(rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_Data,p_DData)
      CALL arrlst_copyArraylist(rarraylist,itable,p_DData)
      
    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Data,p_SData)
      CALL arrlst_copyArraylist(rarraylist,itable,p_SData)
      
    CASE (ST_INT)
      CALL storage_getbase_int(h_Data,p_IData)
      CALL arrlst_copyArraylist(rarraylist,itable,p_IData)
      
    CASE DEFAULT
      PRINT *, 't_arraylist_copyfrom: Unsupported data format!'
      STOP
    END SELECT
  END SUBROUTINE t_arraylist_copyfrom

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyfrom_table(rarraylist,h_Data,h_Table)

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
    IF (h_Table == ST_NOHANDLE) THEN
      CALL storage_new('t_arraylist_copyfrom_table','Table',&
          rarraylist%NTABLE+1,ST_INT,h_Table,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Table,isize)
      IF (isize < rarraylist%NTABLE+1) THEN
        CALL storage_realloc('t_arraylist_copyfrom_table',&
            rarraylist%NTABLE+1,h_Table,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF
    CALL storage_getbase_int(h_Table,p_Table)
    
    IF (h_Data  == ST_NOHANDLE) THEN
      CALL storage_new('t_arraylist_copyfrom_table','Data',&
          rarraylist%NA,rarraylist%carraylistFormat,h_Data,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Data,isize)
      IF (isize < rarraylist%NA) THEN
        CALL storage_realloc('t_arraylist_copyfrom_table',&
            rarraylist%NA,h_Data,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF

    ! What kind of array are we?
    SELECT CASE(rarraylist%carraylistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_Data,p_DData)
      CALL arrlst_copyArraylistTable(rarraylist,p_DData,p_Table)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Data,p_SData)
      CALL arrlst_copyArraylistTable(rarraylist,p_SData,p_Table)

    CASE (ST_INT)
      CALL storage_getbase_int(h_Data,p_IData)
      CALL arrlst_copyArraylistTable(rarraylist,p_IData,p_Table)
      
    CASE DEFAULT
      PRINT *, 't_arraylist_copyfrom_table: Unsupported data format!'
      STOP
    END SELECT
  END SUBROUTINE t_arraylist_copyfrom_table
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyfromDble(rarraylist,itable,p_DData,ndata)

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
      PRINT *, "t_arraylist_copyfromDble: Invalid table number"
      STOP
    END IF
    
    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_copyfromDble: Unsupported data format!"
      STOP
    END IF

    icount = 0
    ipos = rarraylist%Ktable(LHEAD,itable)
    DO
      icount = icount+1
      p_DData(icount) = rarraylist%DData(ipos)
      IF (ipos == rarraylist%Ktable(LTAIL,itable)) EXIT
      ipos = rarraylist%Knext(ipos)
    END DO

    IF (PRESENT(ndata)) ndata=icount
  END SUBROUTINE t_arraylist_copyfromDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyfromDble_table(rarraylist,p_DData,p_Table)

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
    IF (ntable /= rarraylist%NTABLE) THEN
      PRINT *, "t_arraylist_copyfromDble_table: Invalid dimension of table array!"
      STOP
    END IF
    
    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_copyfromDble_table: Unsupported data format!"
      STOP
    END IF

    icount=1
    DO itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%Ktable(LHEAD,itable)
      DO WHILE (ipos /= ANULL)
        p_DData(icount) = rarraylist%DData(ipos)
        icount = icount+1
        IF (ipos == rarraylist%Ktable(LTAIL,itable)) EXIT
        ipos = rarraylist%Knext(ipos)
      END DO
    END DO
    p_Table(ntable+1)=icount+1
  END SUBROUTINE t_arraylist_copyfromDble_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyfromSngl(rarraylist,itable,p_SData,ndata)

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
      PRINT *, "t_arraylist_copyfromSngl: Invalid table number"
      STOP
    END IF

    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "_arraylist_copyfromSngl: Unsupported data format!"
      STOP
    END IF

    icount = 0
    ipos = rarraylist%Ktable(LHEAD,itable)
    DO
      icount = icount+1
      p_SData(icount) = rarraylist%SData(ipos)
      IF (ipos == rarraylist%Ktable(LTAIL,itable)) EXIT
      ipos = rarraylist%Knext(ipos)
    END DO

    IF (PRESENT(ndata)) ndata=icount
  END SUBROUTINE t_arraylist_copyfromSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyfromSngl_table(rarraylist,p_SData,p_Table)

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
    IF (ntable /= rarraylist%NTABLE) THEN
      PRINT *, "t_arraylist_copyfromSngl_table: Invalid dimension of table array!"
      STOP
    END IF
    
    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "t_arraylist_copyfromSngl_table: Unsupported data format!"
      STOP
    END IF

    icount=1
    DO itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%Ktable(LHEAD,itable)
      DO WHILE(ipos /= ANULL)
        p_SData(icount) = rarraylist%SData(ipos)
        icount = icount+1
        IF (ipos == rarraylist%Ktable(LTAIL,itable)) EXIT
        ipos = rarraylist%Knext(ipos)
      END DO
    END DO
    p_Table(ntable+1)=icount+1
  END SUBROUTINE t_arraylist_copyfromSngl_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyfromInt(rarraylist,itable,p_IData,ndata)

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
      PRINT *, "t_arraylist_copyfromInt: Invalid table number"
      STOP
    END IF

    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_copyfromInt: Unsupported data format!"
      STOP
    END IF

    icount = 0
    ipos = rarraylist%Ktable(LHEAD,itable)
    DO
      icount = icount+1
      p_IData(icount) = rarraylist%IData(ipos)
      IF (ipos == rarraylist%Ktable(LTAIL,itable)) EXIT
      ipos = rarraylist%Knext(ipos)
    END DO

    IF (PRESENT(ndata)) ndata=icount
  END SUBROUTINE t_arraylist_copyfromInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyfromInt_table(rarraylist,p_IData,p_Table)

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
    IF (ntable /= rarraylist%NTABLE) THEN
      PRINT *, "t_arraylist_copyfromInt_table: Invalid dimension of table array!"
      STOP
    END IF
    
    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_copyfromInt_table: Unsupported data format!"
      STOP
    END IF

    icount=1
    DO itable=1,ntable
      p_Table(itable) = icount
      
      ipos = rarraylist%Ktable(LHEAD,itable)
      DO WHILE(ipos /= ANULL)
        p_IData(icount) = rarraylist%IData(ipos)
        icount=icount+1
        IF (ipos == rarraylist%Ktable(LTAIL,itable)) EXIT
        ipos = rarraylist%Knext(ipos)
      END DO
    END DO
    p_Table(ntable+1)=icount
  END SUBROUTINE t_arraylist_copyfromInt_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyto(h_DataSrc,itable,rarraylist)

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
      CALL arrlst_copyArraylist(p_DData,itable,rarraylist)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_DataSrc,p_SData)
      CALL arrlst_copyArraylist(p_SData,itable,rarraylist)

    CASE (ST_INT)
      CALL storage_getbase_int(h_DataSrc,p_IData)
      CALL arrlst_copyArraylist(p_IData,itable,rarraylist)
      stop

    CASE DEFAULT
      PRINT *, "t_arraylist_copy: Unsupported data format!"
      STOP
    END SELECT
  END SUBROUTINE t_arraylist_copyto

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copyto_table(h_DataSrc,rarraylist,h_Table)

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
      CALL arrlst_copyArraylistTable(p_DData,rarraylist,p_Table)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_DataSrc,p_SData)
      CALL arrlst_copyArraylistTable(p_SData,rarraylist,p_Table)

    CASE (ST_INT)
      CALL storage_getbase_int(h_DataSrc,p_IData)
      CALL arrlst_copyArraylistTable(p_IData,rarraylist,p_Table)

    CASE DEFAULT
      PRINT *, "t_arraylist_copy_table: Unsupported data format!"
      STOP
    END SELECT
  END SUBROUTINE t_arraylist_copyto_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copytoDble(p_DDataSrc,itable,rarraylist)

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

    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_copytoDble: Unsupported data format!"
      STOP
    END IF
    
    DO ipos=1,SIZE(p_DDataSrc)
      CALL append(rarraylist,itable,p_DDataSrc(ipos),kpos)
    END DO
  END SUBROUTINE t_arraylist_copytoDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copytoDble_table(p_DDataSrc,rarraylist,p_Table)

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

    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_copytoDble_table: Unsupported data format!"
      STOP
    END IF
    
    ntable=SIZE(p_Table)-1
    DO itable=1,ntable
      DO ipos=p_Table(itable),p_Table(itable+1)-1
        CALL append(rarraylist,itable,p_DDataSrc(ipos),kpos)
      END DO
    END DO
  END SUBROUTINE t_arraylist_copytoDble_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copytoSngl(p_SDataSrc,itable,rarraylist)

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

    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "t_arraylist_copytoSngl: Unsupported data format!"
      STOP
    END IF
    
    DO ipos=1,SIZE(p_SDataSrc)
      CALL append(rarraylist,itable,p_SDataSrc(ipos),kpos)
    END DO
  END SUBROUTINE t_arraylist_copytoSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copytoSngl_table(p_SDataSrc,rarraylist,p_Table)

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

    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "t_arraylist_copytoSngl_table: Unsupported data format!"
      STOP
    END IF

    ntable=SIZE(p_Table)-1
    DO itable=1,ntable
      DO ipos=p_Table(itable),p_Table(itable+1)-1
        CALL append(rarraylist,itable,p_SDataSrc(ipos),kpos)
      END DO
    END DO
  END SUBROUTINE t_arraylist_copytoSngl_table

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copytoInt(p_IDataSrc,itable,rarraylist)

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

    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_copytoInt: Unsupported data format!"
      STOP
    END IF
    
    DO ipos=1,SIZE(p_IDataSrc)
      CALL append(rarraylist,itable,p_IDataSrc(ipos),kpos)
    END DO
  END SUBROUTINE t_arraylist_copytoInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_copytoInt_table(p_IDataSrc,rarraylist,p_Table)

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

    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_copytoInt_table: Unsupported data format!"
      STOP
    END IF
    
    ntable=SIZE(p_Table)-1
    DO itable=1,ntable
      DO ipos=p_Table(itable),p_Table(itable+1)-1
        CALL append(rarraylist,itable,p_IDataSrc(ipos),kpos)
      END DO
    END DO
  END SUBROUTINE t_arraylist_copytoInt_table
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_arraylist_swap(rarraylist,itable,jtable)

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
    INTEGER(PREC_ARRAYLISTIDX) :: ihead,itail,iitem
    
    ! Swap
    ihead = rarraylist%Ktable(LHEAD,itable)
    itail = rarraylist%Ktable(LTAIL,itable)
    iitem = rarraylist%Ktable(LITEM,itable)
    
    rarraylist%Ktable(LHEAD,itable) = rarraylist%Ktable(LHEAD,jtable)
    rarraylist%Ktable(LTAIL,itable) = rarraylist%Ktable(LTAIL,jtable)
    rarraylist%Ktable(LITEM,itable) = rarraylist%Ktable(LITEM,jtable)

    rarraylist%Ktable(LHEAD,jtable) = ihead
    rarraylist%Ktable(LTAIL,jtable) = itail
    rarraylist%Ktable(LITEM,jtable) = iitem
  END SUBROUTINE t_arraylist_swap

  ! ***************************************************************************
  
!<function>
  
  PURE FUNCTION t_arraylist_first(rarraylist,itable) RESULT(ipos)

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
    
    ipos=rarraylist%Knext(rarraylist%Ktable(LHEAD,itable))
  END FUNCTION t_arraylist_first

  ! ***************************************************************************
  
!<function>
  
  PURE FUNCTION t_arraylist_last(rarraylist,itable) RESULT(ipos)

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
    
    ipos=rarraylist%Knext(rarraylist%Ktable(LTAIL,itable))
  END FUNCTION t_arraylist_last

  ! ***************************************************************************

!<function>

  FUNCTION t_arraylist_next(rarraylist,itable,breset) RESULT(ipos)

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
      ipos = rarraylist%Ktable(LHEAD,itable)
      rarraylist%Ktable(LITEM,itable) = ipos
      RETURN
    END IF

    ! Get next item and increase item pointer
    ipos = rarraylist%Knext(rarraylist%Ktable(LITEM,itable))
    rarraylist%Ktable(LITEM,itable) = ipos
  END FUNCTION t_arraylist_next

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_prependDble(rarraylist,itable,da,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_prependDble: Unsupported data format!"
      STOP
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_prependDble: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)
    
    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
   
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%DData(ipos)          = da
    ELSE
      rarraylist%Knext(ipos)          = rarraylist%Ktable(LHEAD,itable)
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%DData(ipos)          = da
    END IF
  END SUBROUTINE t_arraylist_prependDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_prependSngl(rarraylist,itable,sa,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "t_arraylist_prependSngl: Unsupported data format!"
      STOP
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_prependSngl: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%SData(ipos)          = sa
    ELSE
      rarraylist%Knext(ipos)          = rarraylist%Ktable(LHEAD,itable)
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%SData(ipos)          = sa
    END IF
  END SUBROUTINE t_arraylist_prependSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_prependInt(rarraylist,itable,ia,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_prependInt: Unsupported data format!"
      STOP
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_prependInt: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%IData(ipos)          = ia
    ELSE
      rarraylist%Knext(ipos)          = rarraylist%Ktable(LHEAD,itable)
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%IData(ipos)          = ia
    END IF
  END SUBROUTINE t_arraylist_prependInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_appendDble(rarraylist,itable,da,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_appendDble: Unsupported data format!"
      STOP
    END IF

    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_appendDble: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)
    
    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%DData(ipos)          = da
    ELSE
      rarraylist%Knext(rarraylist%Ktable(LTAIL,itable)) = ipos
      rarraylist%Ktable(LTAIL,itable)                   = ipos
      rarraylist%Knext(ipos)                            = ANULL
      rarraylist%DData(ipos)                            = da
    END IF
  END SUBROUTINE t_arraylist_appendDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_appendSngl(rarraylist,itable,sa,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "t_arraylist_appendSngl: Unsupported data format!"
      STOP
    END IF

    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_appendSngl: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%SData(ipos)          = sa
    ELSE
      rarraylist%Knext(rarraylist%Ktable(LTAIL,itable)) = ipos
      rarraylist%Ktable(LTAIL,itable)                   = ipos
      rarraylist%Knext(ipos)                            = ANULL
      rarraylist%SData(ipos)                            = sa
    END IF
  END SUBROUTINE t_arraylist_appendSngl
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_appendInt(rarraylist,itable,ia,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_appendInt: Unsupported data format!"
      STOP
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_appendInt: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%IData(ipos)          = ia
    ELSE
      rarraylist%Knext(rarraylist%Ktable(LTAIL,itable)) = ipos
      rarraylist%Ktable(LTAIL,itable)                   = ipos
      rarraylist%Knext(ipos)                            = ANULL
      rarraylist%IData(ipos)                            = ia
    END IF
  END SUBROUTINE t_arraylist_appendInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_insertDble(rarraylist,itable,da,ipred,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_insertDble: Unsupported data format!"
      STOP
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_insertDble: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%DData(ipos)          = da
    ELSEIF (ipred == rarraylist%Ktable(LTAIL,itable)) THEN
      rarraylist%Knext(ipred)         = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%DData(ipos)          = da
    ELSEIF (ipred < ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Knext(ipos)          = -ipred
      rarraylist%DData(ipos)          = da
    ELSE
      rarraylist%Knext(ipos)          = rarraylist%Knext(ipred)
      rarraylist%Knext(ipred)         = ipos
      rarraylist%DData(ipos)          = da
    END IF
  END SUBROUTINE t_arraylist_insertDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_insertSngl(rarraylist,itable,sa,ipred,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "t_arraylist_insertSngl: Unsupported data format!"
      STOP
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_insertSngl: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%SData(ipos)          = sa
    ELSEIF (ipred == rarraylist%Ktable(LTAIL,itable)) THEN
      rarraylist%Knext(ipred)         = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%SData(ipos)          = sa
    ELSEIF (ipred < ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Knext(ipos)          = -ipred
      rarraylist%SData(ipos)          = sa
    ELSE
      rarraylist%Knext(ipos)          = rarraylist%Knext(ipred)
      rarraylist%Knext(ipred)         = ipos
      rarraylist%SData(ipos)          = sa
    END IF
  END SUBROUTINE t_arraylist_insertSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_insertInt(rarraylist,itable,ia,ipred,ipos)

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
    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_insertInt: Unsupported data format!"
      STOP
    END IF
    
    ! Check if tables need to be created
    IF (itable < 1) THEN
      PRINT *, "t_arraylist_insertInt: Invalid table number!"
      STOP
    END IF
    IF (rarraylist%NTABLE < itable) CALL create_table(rarraylist,itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos = rarraylist%Knext(LFREE)
    IF (ABS(ipos) > rarraylist%NNA) THEN
      CALL resize(rarraylist,CEILING(rarraylist%dfactor*rarraylist%NNA))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rarraylist%Knext(LFREE) = ipos+1
    ELSE
      ipos = ABS(ipos)
      rarraylist%Knext(LFREE) = rarraylist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    IF (rarraylist%Ktable(LHEAD,itable) == ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%IData(ipos)          = ia
    ELSEIF (ipred == rarraylist%Ktable(LTAIL,itable)) THEN
      rarraylist%Knext(ipred)         = ipos
      rarraylist%Ktable(LTAIL,itable) = ipos
      rarraylist%Knext(ipos)          = ANULL
      rarraylist%IData(ipos)          = ia
    ELSEIF (ipred < ANULL) THEN
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Knext(ipos)          = -ipred
      rarraylist%IData(ipos)          = ia
    ELSE
      rarraylist%Knext(ipos)          = rarraylist%Knext(ipred)
      rarraylist%Knext(ipred)         = ipos
      rarraylist%IData(ipos)          = ia
    END IF
  END SUBROUTINE t_arraylist_insertInt
  
  ! ***************************************************************************
  
!<function>

  FUNCTION t_arraylist_deleteDble(rarraylist,itable,da) RESULT(f)

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
    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_deleteDble: Unsupported data format!"
      STOP
    END IF

    ! Search for data
    f=search(rarraylist,itable,da,ipred)
    IF (f == ARRAYLIST_NOT_FOUND) RETURN
    
    ! Delete data
    rarraylist%NA = rarraylist%NA-1

    ! Are we first entry in list?
    IF (ipred < 0) THEN

      ! Get position
      ipred = -ipred
      ipos  = rarraylist%Knext(ipred)
      
      ! Update free position
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Knext(ipred)         = rarraylist%Knext(LFREE)
      rarraylist%Knext(LFREE)         = -ipred

    ELSE

      ! Get position
      ipos = rarraylist%Knext(ipred)
      IF (rarraylist%Knext(ipred) == rarraylist%Ktable(LTAIL,itable))&
          rarraylist%Ktable(LTAIL,itable)=ipred
      
      ! Update free position
      rarraylist%Knext(ipred) = rarraylist%Knext(ipos)
      rarraylist%Knext(ipos)  = rarraylist%Knext(LFREE)
      rarraylist%Knext(LFREE) = -ipos
    END IF
  END FUNCTION t_arraylist_deleteDble
  
  ! ***************************************************************************
  
!<function>

  FUNCTION t_arraylist_deleteSngl(rarraylist,itable,sa) RESULT(f)

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
    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "t_arraylist_deleteSngl: Unsupported data format!"
      STOP
    END IF

    ! Search for data
    f=search(rarraylist,itable,sa,ipred)
    IF (f == ARRAYLIST_NOT_FOUND) RETURN

    ! Delete data
    rarraylist%NA = rarraylist%NA-1
    
    ! Are we first entry in list?
    IF (ipred < 0) THEN
      
      ! Get position
      ipred = -ipred
      ipos  = rarraylist%Knext(ipred)
      
      ! Update free position
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Knext(ipred)         = rarraylist%Knext(LFREE)
      rarraylist%Knext(LFREE)         = -ipred
      
    ELSE
      
      ! Get position
      ipos = rarraylist%Knext(ipred)
      IF (rarraylist%Knext(ipred) == rarraylist%Ktable(LTAIL,itable))&
          rarraylist%Ktable(LTAIL,itable)=ipred
      
      ! Update free position
      rarraylist%Knext(ipred) = rarraylist%Knext(ipos)
      rarraylist%Knext(ipos)  = rarraylist%Knext(LFREE)
      rarraylist%Knext(LFREE) = -ipos
    END IF
  END FUNCTION t_arraylist_deleteSngl

  ! ***************************************************************************
  
!<function>

  FUNCTION t_arraylist_deleteInt(rarraylist,itable,ia) RESULT(f)

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
    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_deleteInt: Unsupported data format!"
      STOP
    END IF

    ! Search for data
    f=search(rarraylist,itable,ia,ipred)
    IF (f == ARRAYLIST_NOT_FOUND) RETURN

    ! Delete data
    rarraylist%NA = rarraylist%NA-1

    ! Are we first entry in list?
    IF (ipred < 0) THEN

      ! Get position
      ipred = -ipred
      ipos  = rarraylist%Knext(ipred)

      ! Update free position
      rarraylist%Ktable(LHEAD,itable) = ipos
      rarraylist%Knext(ipred)         = rarraylist%Knext(LFREE)
      rarraylist%Knext(LFREE)         = -ipred
      
    ELSE

      ! Get position
      ipos = rarraylist%Knext(ipred)

      ! Check if last entry should be deleted
      IF (rarraylist%Knext(ipred) == rarraylist%Ktable(LTAIL,itable))&
          rarraylist%Ktable(LTAIL,itable)=ipred
      
      ! Update free position
      rarraylist%Knext(ipred) = rarraylist%Knext(ipos)
      rarraylist%Knext(ipos)  = rarraylist%Knext(LFREE)
      rarraylist%Knext(LFREE) = -ipos
    END IF
  END FUNCTION t_arraylist_deleteInt

  ! ***************************************************************************
  
!<function>

  FUNCTION t_arraylist_searchDble(rarraylist,itable,da,ipred) RESULT(f)

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
    IF (rarraylist%carraylistFormat /= ST_DOUBLE) THEN
      PRINT *, "t_arraylist_searchDble: Unsupported data format!"
      STOP
    END IF

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) RETURN
    ipred=-rarraylist%Ktable(LHEAD,itable)
    
    ! Check if list is empty
    IF (ipred == ANULL) RETURN

    ! Initialization
    ihead=rarraylist%Ktable(LHEAD,itable)
    itail=rarraylist%Ktable(LTAIL,itable)

    ! What kind of ordering are we
    SELECT CASE(rarraylist%cordering)
    CASE (ARRAYLIST_UNORDERED)

      ! Check first item separately
      IF (rarraylist%DData(ihead) == da) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%DData(inext) == da) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (inext == itail) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_INCREASING)

      ! Check first item separately
      IF (rarraylist%DData(ihead) == da) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%DData(ihead) > da) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%DData(inext) == da) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%DData(inext) > da) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_DECREASING)

      ! Check first item separately
      IF (rarraylist%DData(ihead) == da) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%DData(ihead) < da) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)     
        IF (rarraylist%DData(inext) == da) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%DData(inext) < da) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_CSR7)

      ! Check first item separately
      IF (rarraylist%DData(ihead) == da) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred
      
      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%DData(inext) == da) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%DData(inext) > da) EXIT
        IF (rarraylist%Knext(ipred) == itail) THEN
          ipred=rarraylist%Knext(ipred); EXIT
        END IF
        ipred=rarraylist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION t_arraylist_searchDble
  
  ! ***************************************************************************
  
!<function>

  FUNCTION t_arraylist_searchSngl(rarraylist,itable,sa,ipred) RESULT(f)

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
    IF (rarraylist%carraylistFormat /= ST_SINGLE) THEN
      PRINT *, "t_arraylist_searchSngl: Unsupported data format!"
      STOP
    END IF

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) RETURN
    ipred=-rarraylist%Ktable(LHEAD,itable)

    ! Check if list is empty
    IF (ipred == ANULL) RETURN

    ! Initialization
    ihead=rarraylist%Ktable(LHEAD,itable)
    itail=rarraylist%Ktable(LTAIL,itable)

    ! What kind of ordering are we
    SELECT CASE(rarraylist%cordering)
    CASE (ARRAYLIST_UNORDERED)

      ! Check first item separately
      IF (rarraylist%SData(ihead) == sa) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%SData(inext) == sa) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (inext == itail) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_INCREASING)

      ! Check first item separately
      IF (rarraylist%SData(ihead) == sa) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%SData(ihead) > sa) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%SData(inext) == sa) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%SData(inext) > sa) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_DECREASING)

      ! Check first item separately
      IF (rarraylist%SData(ihead) == sa) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%SData(ihead) < sa) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%SData(inext) == sa) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%SData(inext) < sa) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_CSR7)

      ! Check first item separately
      IF (rarraylist%SData(ihead) == sa) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%SData(inext) == sa) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%SData(inext) > sa) EXIT
        IF (rarraylist%Knext(ipred) == itail) THEN
          ipred=rarraylist%Knext(ipred); EXIT
        END IF
        ipred=rarraylist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION t_arraylist_searchSngl

  ! ***************************************************************************
  
!<function>

  FUNCTION t_arraylist_searchInt(rarraylist,itable,ia,ipred) RESULT(f)

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
    IF (rarraylist%carraylistFormat /= ST_INT) THEN
      PRINT *, "t_arraylist_searchInt: Unsupported data format!"
      STOP
    END IF

    ! Initialization
    f=ARRAYLIST_NOT_FOUND

    ! Check if table exists
    IF (itable < 1 .OR. itable > rarraylist%NTABLE) RETURN
    ipred=-rarraylist%Ktable(LHEAD,itable)
    
    ! Check if list is empty
    IF (ipred == ANULL) RETURN

    ! Initialization
    ihead=rarraylist%Ktable(LHEAD,itable)
    itail=rarraylist%Ktable(LTAIL,itable)
        
    ! What kind of ordering are we
    SELECT CASE(rarraylist%cordering)
    CASE (ARRAYLIST_UNORDERED)

      ! Check first item separately
      IF (rarraylist%IData(ihead) == ia) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%IData(inext) == ia) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF

        IF (inext == itail) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_INCREASING)

      ! Check first item separately
      IF (rarraylist%IData(ihead) == ia) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%IData(ihead) > ia) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%IData(inext) == ia) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%IData(inext) > ia) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_DECREASING)

      ! Check first item separately
      IF (rarraylist%IData(ihead) == ia) THEN
        f=ARRAYLIST_FOUND; RETURN
      ELSEIF(rarraylist%IData(ihead) < ia) THEN
        RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%IData(inext) == ia) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%IData(inext) < ia) EXIT
        ipred=rarraylist%Knext(ipred)
      END DO
      
    CASE (ARRAYLIST_CSR7)

      ! Check first item separately
      IF (rarraylist%IData(ihead) == ia) THEN
        f=ARRAYLIST_FOUND; RETURN
      END IF
      ipred=-ipred

      DO WHILE(ipred.NE.itail)
        inext = rarraylist%Knext(ipred)
        IF (rarraylist%IData(inext) == ia) THEN
          f=ARRAYLIST_FOUND; EXIT
        END IF
        
        IF (rarraylist%IData(inext) > ia) EXIT
        IF (rarraylist%Knext(ipred) == itail) THEN
          ipred=rarraylist%Knext(ipred); EXIT
        END IF
        ipred=rarraylist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION t_arraylist_searchInt

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE t_arraylist_print(rarraylist,itable)

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
      ! Loop over all selected tables
      DO iitable=itable1,itable2
        
        WRITE(*,FMT='(A,I5)') 'Table No.:',iitable
        WRITE(*,FMT='(A)',ADVANCE='NO') ' ->'

        ipos  = rarraylist%Ktable(LHEAD,iitable)
        IF (ipos == ANULL) CYCLE
        itail = rarraylist%Ktable(LTAIL,iitable)
        
        DO
          WRITE(*,FMT='(I5,1X)',ADVANCE='NO') rarraylist%DData(ipos)
          IF (ipos == itail) EXIT
          ipos = rarraylist%Knext(ipos)
        END DO
        PRINT *
      END DO
      
    CASE (ST_SINGLE)
      ! Loop over all selected tables
      DO iitable=itable1,itable2

        WRITE(*,FMT='(A,I5)') 'Table No.:',iitable
        WRITE(*,FMT='(A)',ADVANCE='NO') ' ->'
        
        ipos = rarraylist%Ktable(LHEAD,iitable)
        IF (ipos == ANULL) CYCLE
        itail = rarraylist%Ktable(LTAIL,iitable)
        
        DO
          WRITE(*,FMT='(I5,1X)',ADVANCE='NO') rarraylist%SData(ipos)
          IF (ipos == itail) EXIT
          ipos = rarraylist%Knext(ipos)
        END DO
        PRINT *
      END DO
      
    CASE (ST_INT)
      ! Loop over all selected tables
      DO iitable=itable1,itable2

        WRITE(*,FMT='(A,I5)') 'Table No.:',iitable
        WRITE(*,FMT='(A)',ADVANCE='NO') ' ->'

        ipos = rarraylist%Ktable(LHEAD,iitable)
        IF (ipos == ANULL) CYCLE
        itail = rarraylist%Ktable(LTAIL,iitable)
        
        DO
          WRITE(*,FMT='(I5,1X)',ADVANCE='NO') rarraylist%IData(ipos)
          IF (ipos == itail) EXIT
          ipos = rarraylist%Knext(ipos)
        END DO
        PRINT *
      END DO
      
    CASE DEFAULT
      PRINT *, "t_arraylist_print: Unsupported data type!"
      STOP
    END SELECT
  END SUBROUTINE t_arraylist_print

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_arraylist_info(rarraylist)

!<description>
    ! This subroutine prints information about the arraylist
!</description>

!<input>
    ! arraylist
    TYPE(t_arraylist), INTENT(IN) :: rarraylist
!</input>
!</subroutine>

    WRITE(*,FMT=*) ' Arraylist:'
    WRITE(*,FMT=*) ' =========='
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Ktable =',rarraylist%h_Ktable
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Knext  =',rarraylist%h_Knext
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Data   =',rarraylist%h_Data
    WRITE(*,*)
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8,3X,A,3X,F5.1,A)') &
        '  NA      =',rarraylist%NA,&
        'NNA      =',rarraylist%NNA,&
        'FILLING  =',100*rarraylist%NA/REAL(rarraylist%NNA,DP),'%'
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8,3X,A,3X,F5.1,A)') &
        '  NTABLE  =',rarraylist%NTABLE,&
        'NNTABLE  =',rarraylist%NNTABLE,&
        'FILLING  =',100*rarraylist%NTABLE/REAL(rarraylist%NNTABLE,DP)
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,2X,F7.1,A,1X,F7.1,A)') &
        '  NRESIZE =',rarraylist%NRESIZE,&
        'TABLE´s  =',100*rarraylist%NNTABLE/REAL(rarraylist%NNTABLE0,DP),&
        '%  NNA´s    =',100*rarraylist%NNA/REAL(rarraylist%NNA0,DP),'%'
    WRITE(*,*)
  END SUBROUTINE t_arraylist_info
END MODULE arraylist
