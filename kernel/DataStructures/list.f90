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

MODULE list
  USE fsystem
  USE storage
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_list
  PUBLIC :: list_createList
  PUBLIC :: list_releaseList
  PUBLIC :: list_resizeList
  PUBLIC :: list_copyList
  PUBLIC :: list_swapList
  PUBLIC :: list_getFirstInList
  PUBLIC :: list_getLastInList
  PUBLIC :: list_getNextInList
  PUBLIC :: list_getPrevInList
  PUBLIC :: list_prependToList
  PUBLIC :: list_appendToList
  PUBLIC :: list_insertIntoList
  PUBLIC :: list_deleteFromList
  PUBLIC :: list_searchInList
  PUBLIC :: list_printList

!<constants>

!<constantblock description="KIND values for list data">

  ! kind value for indices in list
  INTEGER, PARAMETER, PUBLIC :: PREC_LISTIDX = I32

!</constantblock>

!<constantblock description="Global flags for list ordering">

  ! Identifier for unordered list
  INTEGER, PARAMETER, PUBLIC :: LIST_UNORDERED  = 0
  
  ! Identifier for increasingly ordered list
  INTEGER, PARAMETER, PUBLIC :: LIST_INCREASING = 1

  ! Identifier for decreasingly ordered list
  INTEGER, PARAMETER, PUBLIC :: LIST_DECREASING = 2

  ! Identifier for ordered list 
  INTEGER, PARAMETER, PUBLIC :: LIST_CSR7       = 3

!</constantblock>

!<constantblock description="Global flags for list likn-type">

  ! Identifier fir single-linked list
  INTEGER, PARAMETER, PUBLIC :: LIST_SINGLELINKED = 1
  
  ! Identifier for double-linked list
  INTEGER, PARAMETER, PUBLIC :: LIST_DOUBLELINKED = 2

!</constantblock>

!<constantblock description="Global flags for list operations">

  ! Identifier for "not found in list"
  INTEGER, PARAMETER, PUBLIC :: LIST_NOT_FOUND = -1

  ! Identifier for "found in list"
  INTEGER, PARAMETER, PUBLIC :: LIST_FOUND     =  0
  
!</constantblock>

!<constantblock description="Internal tags for list status">
  
  ! Tag for empty list
  INTEGER(PREC_LISTIDX), PARAMETER, PUBLIC :: LNULL =  0
  
  ! Tag for head of list
  INTEGER(PREC_LISTIDX), PARAMETER :: LHEAD = -2

  ! Tag for tail of list
  INTEGER(PREC_LISTIDX), PARAMETER :: LTAIL = -1

  ! Tag for next free item
  INTEGER(PREC_LISTIDX), PARAMETER :: LFREE =  0

!</constantblock>
!</constants>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<types>
!<typeblock>

  TYPE t_list
    ! Format-Tag: Double, Single, Integer
    INTEGER :: clistFormat        = ST_NOHANDLE
    
    ! Type-Tag: Single-linked, Double-linked
    INTEGER :: clinkType          = LIST_UNORDERED

    ! Type of list ordering
    INTEGER :: cordering          = LIST_UNORDERED
    
    ! Position of the last item inserted into the list
    INTEGER :: item

    ! Number of items that are currently stored in the list
    INTEGER(PREC_LISTIDX) :: NA   = 0

    ! Total number of items that can be stored in the list
    INTEGER(PREC_LISTIDX) :: NNA  = 0

    ! Total number of resize operations
    INTEGER :: NRESIZE            = 0

    ! Dimension of the auxiliary Integer values to be stored
    INTEGER :: isizeInt            = 0

    ! Dimension of the auxiliary Double values to be stored
    INTEGER :: isizeDble           = 0

    ! Dimension of the auxiliary Single values to be stored
    INTEGER :: isizeSngl           = 0

    ! Factor by which the list is enlarged if new storage is allocate
    REAL(DP) :: dfactor           = 1.5_DP

    ! Handle to the list key
    INTEGER :: h_Key              = ST_NOHANDLE

    ! Handle to the list next-structure
    INTEGER :: h_Knext            = ST_NOHANDLE

    ! Handle to the list previous-structure
    INTEGER :: h_Kprev            = ST_NOHANDLE

    ! Handle to the list auxiliary Integer data
    INTEGER :: h_IData             = ST_NOHANDLE

    ! Handle to the list auxiliary Double data
    INTEGER :: h_DData             = ST_NOHANDLE

    ! Handle to the list auxiliary Single data
    INTEGER :: h_SData             = ST_NOHANDLE
    
    ! List next-structure
    ! NOTE: This array is introduced to increase performance. It
    ! should not be touched by the user. However, if the handle would
    ! be dereferenced for each operation such as search, delete,
    ! performance would be very poor.
    INTEGER(PREC_LISTIDX), DIMENSION(:), POINTER :: Knext => NULL()

    ! List previous structure
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    INTEGER(PREC_LISTIDX), DIMENSION(:), POINTER :: Kprev => NULL()

    ! List key data (Integer)
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    INTEGER(PREC_LISTIDX), DIMENSION(:), POINTER :: IKey => NULL()

    ! List key data (Double)
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    REAL(DP), DIMENSION(:), POINTER :: DKey => NULL()

    ! List key data (Single)
    ! NOTE: This array is introduced to increase performance (see
    ! above)
    REAL(SP), DIMENSION(:), POINTER :: SKey => NULL()
    
    ! List data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(DP), DIMENSION(:,:), POINTER :: DData => NULL()

    ! List data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(SP), DIMENSION(:,:), POINTER :: SData => NULL()

    ! List data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_LISTIDX), DIMENSION(:,:), POINTER :: IData => NULL()
  END TYPE t_list
  
!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  INTERFACE list_createList
    MODULE PROCEDURE t_list_create
  END INTERFACE
  
  INTERFACE list_releaseList
    MODULE PROCEDURE t_list_release
  END INTERFACE
  
  INTERFACE list_resizeList
    MODULE PROCEDURE t_list_resize
  END INTERFACE
  INTERFACE resize   ! internal use
    MODULE PROCEDURE t_list_resize
  END INTERFACE
  
  INTERFACE list_copyList
    MODULE PROCEDURE t_list_copyto
    MODULE PROCEDURE t_list_copytoDble
    MODULE PROCEDURE t_list_copytoSngl
    MODULE PROCEDURE t_list_copytoInt
    MODULE PROCEDURE t_list_copyfrom
    MODULE PROCEDURE t_list_copyfromDble
    MODULE PROCEDURE t_list_copyfromSngl
    MODULE PROCEDURE t_list_copyfromInt
  END INTERFACE
  
  INTERFACE list_swapList
    MODULE PROCEDURE t_list_swap
  END INTERFACE
  
  INTERFACE list_getFirstInList
    MODULE PROCEDURE t_list_first
  END INTERFACE
  
  INTERFACE list_getLastInList
    MODULE PROCEDURE t_list_last
  END INTERFACE
  
  INTERFACE list_getNextInList
    MODULE PROCEDURE t_list_next
  END INTERFACE

  INTERFACE list_getPrevInList
    MODULE PROCEDURE t_list_prev
  END INTERFACE
  
  INTERFACE list_prependToList
    MODULE PROCEDURE t_list_prependDble
    MODULE PROCEDURE t_list_prependSngl
    MODULE PROCEDURE t_list_prependInt
  END INTERFACE
  INTERFACE prepend   ! internal use
    MODULE PROCEDURE t_list_prependDble
    MODULE PROCEDURE t_list_prependSngl
    MODULE PROCEDURE t_list_prependInt
  END INTERFACE
  
  INTERFACE list_appendToList
    MODULE PROCEDURE t_list_appendDble
    MODULE PROCEDURE t_list_appendSngl
    MODULE PROCEDURE t_list_appendInt
  END INTERFACE
  INTERFACE append   ! internal use
    MODULE PROCEDURE t_list_appendDble
    MODULE PROCEDURE t_list_appendSngl
    MODULE PROCEDURE t_list_appendInt
  END INTERFACE
  
  INTERFACE list_insertIntoList
    MODULE PROCEDURE t_list_insertDble
    MODULE PROCEDURE t_list_insertSngl
    MODULE PROCEDURE t_list_insertInt
  END INTERFACE
  INTERFACE insert   ! internal use
    MODULE PROCEDURE t_list_insertDble
    MODULE PROCEDURE t_list_insertSngl
    MODULE PROCEDURE t_list_insertInt
  END INTERFACE

  INTERFACE list_deleteFromList
    MODULE PROCEDURE t_list_deleteDble
    MODULE PROCEDURE t_list_deleteSngl
    MODULE PROCEDURE t_list_deleteInt
  END INTERFACE
  INTERFACE delete   ! internal use
    MODULE PROCEDURE t_list_deleteDble
    MODULE PROCEDURE t_list_deleteSngl
    MODULE PROCEDURE t_list_deleteInt
  END INTERFACE
  
  INTERFACE list_searchInList
    MODULE PROCEDURE t_list_searchDble
    MODULE PROCEDURE t_list_searchSngl
    MODULE PROCEDURE t_list_searchInt
  END INTERFACE
  INTERFACE search   ! internal use
    MODULE PROCEDURE t_list_searchDble
    MODULE PROCEDURE t_list_searchSngl
    MODULE PROCEDURE t_list_searchInt
  END INTERFACE
  
  INTERFACE list_printList
    MODULE PROCEDURE t_list_print
  END INTERFACE
  
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

CONTAINS
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE t_list_create(rlist,nna,clistFormat,&
      isizeInt,isizeDble,isizeSngl,cordering,dfactor,clinkType)

!<description>
    ! This subroutine creates a new list
!</description>

!<input>
    ! Total number of items that can be stored in list
    INTEGER(PREC_LISTIDX), INTENT(IN) :: nna

    ! Format-tag. Type of list format (Double,Single,Integer)
    INTEGER, INTENT(IN) :: clistFormat

    ! Dimension of the auxiliary Integer values to be stored
    INTEGER, INTENT(IN) :: isizeInt

    ! Dimension of the auxiliary Double values to be stored
    INTEGER, INTENT(IN) :: isizeDble

    ! Dimension of the auxiliary Single values to be stored
    INTEGER, INTENT(IN) :: isizeSngl

    ! OPTIONAL: Format-tag. Type of list ordering 
    INTEGER, INTENT(IN), OPTIONAL :: cordering

    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    REAL(DP), INTENT(IN), OPTIONAL :: dfactor

    ! OPTIONAL: Type of linking (single/double). If not specified
    ! then a single-linked list is generated
    INTEGER, INTENT(IN), OPTIONAL :: clinkType
!</input>

!<output>
    ! list
    TYPE(t_list), INTENT(OUT) :: rlist
!</output>
!</subroutine>
    
    ! Set factor
    IF (PRESENT(dfactor)) THEN
      IF (dfactor > 1_DP) rlist%dfactor=dfactor
    END IF

    ! Set list format
    rlist%clistFormat=clistFormat

    ! Set ordering
    IF (PRESENT(cordering)) THEN
      rlist%cordering=cordering
    END IF

    ! Initialize list
    rlist%isizeInt  = isizeInt
    rlist%isizeSngl = isizeSngl
    rlist%isizeDble = isizeDble
    rlist%item      = LHEAD
    rlist%na        = 0
    rlist%nna       = nna
    IF (PRESENT(clinkType)) THEN
      rlist%clinkType=clinkType
    ELSE
      rlist%clinkType=LIST_SINGLELINKED
    END IF
    
    ! Allocate memory and associate pointers
    CALL storage_new('t_list_create','Knext',-2,nna,ST_INT,&
        rlist%h_Knext,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(rlist%h_Knext,rlist%Knext)
    
    ! Double-linked list?    
    IF (rlist%clinkType == LIST_DOUBLELINKED) THEN
      CALL storage_new('t_list_create','Kprev',nna,ST_INT,&
          rlist%h_Kprev,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(rlist%h_Kprev,rlist%Kprev)
    END IF

    ! Allocate memory for Key
    SELECT CASE(rlist%clistFormat)
    CASE (ST_DOUBLE)
      CALL storage_new('t_list_create','Key',nna,ST_DOUBLE,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double(rlist%h_Key,rlist%DKey)

    CASE (ST_SINGLE)
      CALL storage_new('t_list_create','Key',nna,ST_SINGLE,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_single(rlist%h_Key,rlist%SKey)

    CASE (ST_INT)
      CALL storage_new('t_list_create','Key',nna,ST_INT,&
          rlist%h_Key,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int(rlist%h_Key,rlist%IKey)

    CASE DEFAULT
      PRINT *, 't_list_create: Unsupported data format!'
      CALL sys_halt()
    END SELECT
    
    ! Initialize list structure
    rlist%Knext(LFREE) = 1
    rlist%Knext(LHEAD) = LNULL
    rlist%Knext(LTAIL) = LNULL

    ! Allocate memory for auxiliary data
    IF (isizeDble > 0) THEN
      CALL storage_new('t_list_create','DData',(/isizeDble,nna/),&
          ST_DOUBLE,rlist%h_DData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double2D(rlist%h_DData,rlist%DData)
    END IF

     IF (isizeSngl > 0) THEN
      CALL storage_new('t_list_create','SData',(/isizeSngl,nna/),&
          ST_SINGLE,rlist%h_SData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_single2D(rlist%h_SData,rlist%SData)
    END IF

    IF (isizeInt > 0) THEN
      CALL storage_new('t_list_create','IData',(/isizeInt,nna/),&
          ST_INT,rlist%h_IData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int2D(rlist%h_IData,rlist%IData)
    END IF
  END SUBROUTINE t_list_create
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_list_release(rlist)

!<description>
    ! This subroutine releases an existing list
!</description>

!<inputoutput>
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>
!</subroutine>

    ! Release memory
    IF (rlist%h_Key /= ST_NOHANDLE)   CALL storage_free(rlist%h_Key)
    IF (rlist%h_Knext /= ST_NOHANDLE) CALL storage_free(rlist%h_Knext)
    IF (rlist%h_Kprev /= ST_NOHANDLE) CALL storage_free(rlist%h_Kprev)

    IF (rlist%h_DData /= ST_NOHANDLE) CALL storage_free(rlist%h_DDATA)
    IF (rlist%h_SData /= ST_NOHANDLE) CALL storage_free(rlist%h_SDATA)
    IF (rlist%h_IData /= ST_NOHANDLE) CALL storage_free(rlist%h_IDATA)

    NULLIFY(rlist%Knext,rlist%Kprev,rlist%DKey,rlist%SKey,rlist%IKey)
    NULLIFY(rlist%DData,rlist%SData,rlist%IData)

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
  END SUBROUTINE t_list_release

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_list_resize(rlist,nna)

!<description>
    ! This subroutine reallocates memory for an existing list
!</description>

!<input>
    ! New number of total items that can be stored in the list
    INTEGER(PREC_LISTIDX), INTENT(IN) :: nna
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! Set new size
    rlist%nna=nna

    ! Reallocate structures
    CALL storage_realloc('t_list_resize',-2,nna,rlist%h_Knext,&
        ST_NEWBLOCK_NOINIT,.TRUE.)
    CALL storage_getbase_int(rlist%h_Knext,rlist%Knext)

    IF (rlist%clinkType == LIST_DOUBLELINKED) THEN
      CALL storage_realloc('t_list_resize',nna,rlist%h_Kprev,&
          ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_int(rlist%h_Kprev,rlist%Kprev)
    END IF

    ! Reallocate Key
    CALL storage_realloc('t_list_resize',nna,rlist%h_Key,&
        ST_NEWBLOCK_NOINIT,.TRUE.)
    SELECT CASE(rlist%clistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(rlist%h_Key,rlist%DKey)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(rlist%h_Key,rlist%SKey)

    CASE (ST_INT)
      CALL storage_getbase_int(rlist%h_Key,rlist%IKey)

    CASE DEFAULT
      PRINT *, 't_list_resize: Unsupported data format!'
      CALL sys_halt()
    END SELECT

    ! Reallocate auxiliary data
    IF (rlist%isizeDble > 0) THEN
      CALL storage_realloc('t_list_resize',nna,rlist%h_DData,&
          ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_double2D(rlist%h_DData,rlist%DData)
    END IF

    IF (rlist%isizeSngl > 0) THEN
      CALL storage_realloc('t_list_resize',nna,rlist%h_SData,&
          ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_single2D(rlist%h_SData,rlist%SData)
    END IF

    IF (rlist%isizeInt > 0) THEN
      CALL storage_realloc('t_list_resize',nna,rlist%h_IData,&
          ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_int2D(rlist%h_IData,rlist%IData)
    END IF
  END SUBROUTINE t_list_resize

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_copyfrom(rlist,h_Key)

!<description>
    ! This subroutine copies the content of the list to the given
    ! handle
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist
!</input>

!<inputoutput>
    ! handle to the data
    INTEGER, INTENT(INOUT) :: h_Key
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_DKey
    REAL(SP), DIMENSION(:), POINTER :: p_SKey
    INTEGER,  DIMENSION(:), POINTER :: p_IKey
    INTEGER(PREC_LISTIDX) :: ipos,jpos
    
    ! Transform the content of the list to h_Key
    IF (h_Key /= ST_NOHANDLE) CALL storage_free(h_Key)
    CALL storage_new('t_list_copy','Key',rlist%na,&
        rlist%clistFormat,h_Key,ST_NEWBLOCK_NOINIT)
    SELECT CASE(rlist%clistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_Key,p_DKey)
      ipos = rlist%Knext(LHEAD)
      jpos = 0
      DO
        jpos = jpos+1
        p_DKey(jpos) = rlist%DKey(ipos)
        IF (ipos == rlist%Knext(LTAIL)) EXIT
        ipos = rlist%Knext(ipos)
      END DO
      
    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Key,p_SKey)
      ipos = rlist%Knext(LHEAD)
      jpos = 0
      DO
        jpos = jpos+1
        p_SKey(jpos) = rlist%SKey(ipos)
        IF (ipos == rlist%Knext(LTAIL)) EXIT
        ipos = rlist%Knext(ipos)
      END DO
      
    CASE (ST_INT)
      CALL storage_getbase_int(h_Key,p_IKey)
      ipos = rlist%Knext(LHEAD)
      jpos = 0
      DO
        jpos = jpos+1
        p_IKey(jpos) = rlist%IKey(ipos)
        IF (ipos == rlist%Knext(LTAIL)) EXIT
        ipos = rlist%Knext(ipos)
      END DO
      
    CASE DEFAULT
      PRINT *, 't_list_copy: Unsupported data format!'
      CALL sys_halt()
    END SELECT
  END SUBROUTINE t_list_copyfrom

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_copyfromDble(rlist,p_DKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! double array
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist
!</input>

!<inputoutput>
    ! double array
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p_DKey
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_LISTIDX) :: ipos,jpos

    IF (rlist%clistFormat /= ST_DOUBLE) THEN
      PRINT *, 't_list_copyfromDble: Unsupported data format!'
      CALL sys_halt()
    END IF

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    DO
      jpos = jpos+1
      p_DKey(jpos) = rlist%DKey(ipos)
      IF (ipos == rlist%Knext(LTAIL)) EXIT
      ipos = rlist%Knext(ipos)
    END DO
  END SUBROUTINE t_list_copyfromDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_copyfromSngl(rlist,p_SKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! single array
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist
!</input>

!<inputoutput>
    ! double array
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: p_SKey
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_LISTIDX) :: ipos,jpos

    IF (rlist%clistFormat /= ST_SINGLE) THEN
      PRINT *, 't_list_copyfromSngl: Unsupported data format!'
      CALL sys_halt()
    END IF

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    DO
      jpos = jpos+1
      p_SKey(jpos) = rlist%SKey(ipos)
      IF (ipos == rlist%Knext(LTAIL)) EXIT
      ipos = rlist%Knext(ipos)
    END DO
  END SUBROUTINE t_list_copyfromSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_copyfromInt(rlist,p_IKey)

!<description>
    ! This subroutine copies the content of the list to the given
    ! integer array
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist
!</input>

!<inputoutput>
    ! double array
    INTEGER, DIMENSION(:), INTENT(INOUT) :: p_IKey
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_LISTIDX) :: ipos,jpos

    IF (rlist%clistFormat /= ST_INT) THEN
      PRINT *, 't_list_copyfromInt: Unsupported data format!'
      CALL sys_halt()
    END IF

    ipos = rlist%Knext(LHEAD)
    jpos = 0
    DO
      jpos = jpos+1
      p_IKey(jpos) = rlist%IKey(ipos)
      IF (ipos == rlist%Knext(LTAIL)) EXIT
      ipos = rlist%Knext(ipos)
    END DO
  END SUBROUTINE t_list_copyfromInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_copyto(h_KeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given handle to the list
!</description>

!<input>
    ! handle to the data
    INTEGER, INTENT(IN) :: h_KeySrc
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_DKey
    REAL(SP), DIMENSION(:), POINTER :: p_SKey
    INTEGER,  DIMENSION(:), POINTER :: p_IKey
    INTEGER(PREC_LISTIDX) :: ipos,kpos
    
    ! Transform the content of h_Data to the list
    SELECT CASE (rlist%clistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_KeySrc,p_DKey)
      DO ipos=1,SIZE(p_DKey)
        CALL append(rlist,p_DKey(ipos),kpos)
      END DO
    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_KeySrc,p_SKey)
      DO ipos=1,SIZE(p_SKey)
        CALL append(rlist,p_DKey(ipos),kpos)
      END DO
    CASE (ST_INT)
      CALL storage_getbase_int(h_KeySrc,p_IKey)
      DO ipos=1,SIZE(p_IKey)
        CALL append(rlist,p_IKey(ipos),kpos)
      END DO
    CASE DEFAULT
      PRINT *, 't_list_copy: Unsupported data format!'
      CALL sys_halt()
    END SELECT
  END SUBROUTINE t_list_copyto

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_copytoDble(p_DKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given double array to the list
!</description>

!<input>
    ! handle to the data
    REAL(DP), DIMENSION(:), INTENT(IN) :: p_DKeySrc
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_LISTIDX) :: ipos,kpos

    IF (rlist%clistFormat /= ST_DOUBLE) THEN
      PRINT *, 't_list_copytoDble: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    DO ipos=1,SIZE(p_DKeySrc)
      CALL append(rlist,p_DKeySrc(ipos),kpos)
    END DO
  END SUBROUTINE t_list_copytoDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_copytoSngl(p_SKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given single array to the list
!</description>

!<input>
    ! handle to the data
    REAL(SP), DIMENSION(:), INTENT(IN) :: p_SKeySrc
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_LISTIDX) :: ipos,kpos

    IF (rlist%clistFormat /= ST_SINGLE) THEN
      PRINT *, 't_list_copytoSngl: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    DO ipos=1,SIZE(p_SKeySrc)
      CALL append(rlist,p_SKeySrc(ipos),kpos)
    END DO
  END SUBROUTINE t_list_copytoSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_copytoInt(p_IKeySrc,rlist)

!<description>
    ! This subroutine copies the content of the given integer array to the list
!</description>

!<input>
    ! handle to the data
    INTEGER, DIMENSION(:), INTENT(IN) :: p_IKeySrc
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_LISTIDX) :: ipos,kpos

    IF (rlist%clistFormat /= ST_INT) THEN
      PRINT *, 't_list_copytoInt: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    DO ipos=1,SIZE(p_IKeySrc)
      CALL append(rlist,p_IKeySrc(ipos),kpos)
    END DO
  END SUBROUTINE t_list_copytoInt

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_list_swap(rlist1,rlist2)

!<description>
    ! This subroutine swaps two lists
!</description>
    
!<inputoutput>
    ! First list
    TYPE(t_list), INTENT(INOUT) :: rlist1

    ! Second list
    TYPE(t_list), INTENT(INOUT) :: rlist2
!</inputoutput>
!</subroutine>

    ! local variables
    TYPE(t_list) :: rlist
    
    ! Check if both lists are compatible
    IF (rlist1%clistFormat /= rlist2%clistFormat .OR.&
        rlist1%clinkType   /= rlist2%clinkType .OR. &
        rlist1%isizeInt    /= rlist2%isizeInt .OR.&
        rlist1%isizeDble   /= rlist2%isizeDble .OR.&
        rlist1%isizeSngl   /= rlist2%isizeSngl) THEN
      PRINT *, "t_list_swap: Lists are not compatible"
      CALL sys_halt()
    END IF

    ! Swap
    rlist  = rlist1
    rlist1 = rlist2
    rlist2 = rlist

    ! Reassociated pointers of first list
    NULLIFY(rlist1%Knext,rlist1%Kprev,rlist1%DKey,&
        rlist1%SKey,rlist1%IKey)
    CALL storage_getbase_int(rlist1%h_Knext,rlist1%Knext)
    IF (rlist1%h_Kprev /= ST_NOHANDLE) &
        CALL storage_getbase_int(rlist1%h_Kprev,rlist1%Kprev)

    SELECT CASE(rlist1%clistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(rlist1%h_Key,rlist1%DKey)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(rlist1%h_Key,rlist1%SKey)

    CASE (ST_INT)
      CALL storage_getbase_int(rlist1%h_Key,rlist1%IKey)

    CASE DEFAULT
      PRINT *, 't_list_swap: Unsupported data type!'
      CALL sys_halt()
    END SELECT

    IF (rlist1%isizeDble > 0)&
        CALL storage_getbase_double2D(rlist1%h_DData,rlist1%DData)
    IF (rlist1%isizeSngl > 0)&
        CALL storage_getbase_single2D(rlist1%h_SData,rlist1%SData)
    IF (rlist1%isizeInt > 0)&
        CALL storage_getbase_int2D(rlist1%h_IData,rlist1%IData)
    
    ! Reassociated pointers of second list
    NULLIFY(rlist2%Knext,rlist2%Kprev,rlist2%DKey,&
        rlist2%SKey,rlist2%IKey)
    CALL storage_getbase_int(rlist2%h_Knext,rlist2%Knext)
    IF (rlist2%h_Kprev /= ST_NOHANDLE) &
        CALL storage_getbase_int(rlist2%h_Kprev,rlist2%Kprev)

    SELECT CASE(rlist2%clistFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(rlist2%h_Key,rlist2%DKey)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(rlist2%h_Key,rlist2%SKey)

    CASE (ST_INT)
      CALL storage_getbase_int(rlist2%h_Key,rlist2%IKey)

    CASE DEFAULT
      PRINT *, 't_list_swap: Unsupported data type!'
      CALL sys_halt()
    END SELECT

    IF (rlist2%isizeDble > 0)&
        CALL storage_getbase_double2D(rlist2%h_DData,rlist2%DData)
    IF (rlist2%isizeSngl > 0)&
        CALL storage_getbase_single2D(rlist2%h_SData,rlist2%SData)
    IF (rlist2%isizeInt > 0)&
        CALL storage_getbase_int2D(rlist2%h_IData,rlist2%IData)
  END SUBROUTINE t_list_swap

  ! ***************************************************************************
  
!<function>
  
  PURE FUNCTION t_list_first(rlist) RESULT(ipos)

!<description>
    ! This function returns the position of the first list item
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist
!</input>

!<result>
    ! position of first item
    INTEGER(PREC_LISTIDX) :: ipos
!</result>
!</function>
    
    ipos=rlist%Knext(LHEAD)
  END FUNCTION t_list_first

  ! ***************************************************************************
  
!<function>
  
  PURE FUNCTION t_list_last(rlist) RESULT(ipos)

!<description>
    ! This function returns the position of the last list item
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist
!</input>

!<result>
    ! position of last item
    INTEGER(PREC_LISTIDX) :: ipos
!</result>
!</function>
    
    ipos=rlist%Knext(LTAIL)
  END FUNCTION t_list_last

  ! ***************************************************************************

!<function>

  FUNCTION t_list_next(rlist,breset) RESULT(ipos)

!<description>
    ! This function returns the position of the next list item
    ! and resets the list if required. If the last list item is
    ! reached, then LNULL is returned and the list is reset
    ! implicitely.
!</description>

!<input>
    ! Reset list?
    LOGICAL, INTENT(IN) :: breset
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<result>
    ! position of last item
    INTEGER(PREC_LISTIDX) :: ipos
!</result>
!</function>
    
    ! Reset?
    IF (breset) rlist%item=LHEAD
    
    ipos = rlist%Knext(rlist%item)
    rlist%item = ipos
  END FUNCTION t_list_next

   ! ***************************************************************************

!<function>

  FUNCTION t_list_prev(rlist,breset) RESULT(ipos)

!<description>
    ! This function returns the position of the previous list item
    ! and resets the list if required. This operation is only
    ! available for double-linked lists. If the first list item is
    ! reached, then LNULL is returned and the list is reset
    ! implicitely.
!</description>

!<input>
    ! Reset list?
    LOGICAL, INTENT(IN) :: breset
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<result>
    ! position of previous item
    INTEGER(PREC_LISTIDX) :: ipos
!</result>
!</function>

    IF (rlist%clinkType /= LIST_DOUBLELINKED) THEN
      PRINT *, "t_list_prev: This operation is only available for&
          & double-linked lists!"
      CALL sys_halt()
    END IF
    
    ! Reset?
    IF (breset) rlist%item=LTAIL

    IF (rlist%item == LNULL) THEN
      ipos = rlist%Knext(LTAIL)
    ELSE
      ipos = rlist%Kprev(rlist%item)
    END IF
    rlist%item = ipos
  END FUNCTION t_list_prev

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_prependDble(rlist,dkey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends a Double data to the list
!</description>

!<input>
    ! Double key
    REAL(DP), INTENT(IN) :: dkey

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_DOUBLE) THEN
      PRINT *, 't_list_prependDble: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE (rlist%clinkType)
    CASE (LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      END IF

    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSE
        rlist%Kprev(rlist%Knext(LHEAD)) = ipos
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
        rlist%Kprev(ipos)  = LNULL
      END IF

    CASE DEFAULT
      PRINT *, "t_list_prependDble: Invalid link type!"
      CALL sys_halt()
    END SELECT

    ! Double key
    rlist%Dkey(ipos) = dkey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_prependDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_prependSngl(rlist,skey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends a Single data to the list
!</description>

!<input>
    ! Single key
    REAL(SP), INTENT(IN) :: skey

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_SINGLE) THEN
      PRINT *, 't_list_prependSngl: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE(rlist%clinkType)
    CASE (LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      END IF

    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSE
        rlist%Kprev(rlist%Knext(LHEAD)) = ipos
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
        rlist%Kprev(ipos)  = LNULL
      END IF
    CASE DEFAULT
      PRINT *, "t_list_prependSngl: Invalid linktype!"
      CALL sys_halt()
    END SELECT

    ! Single key
    rlist%Skey(ipos) = skey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_prependSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_prependInt(rlist,ikey,ipos,DData,SData,IData)

!<description>
    ! This subroutine prepends an Integer data to the list
!</description>

!<input>
    ! Integer key
    INTEGER, INTENT(IN) :: ikey

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_INT) THEN
      PRINT *, 't_list_prependInt: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE (rlist%clinkType)
    CASE (LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
      END IF

    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSE
        rlist%Kprev(rlist%Knext(LHEAD)) = ipos
        rlist%Knext(ipos)  = rlist%Knext(LHEAD)
        rlist%Knext(LHEAD) = ipos
        rlist%Kprev(ipos)  = LNULL
      END IF
      
    CASE DEFAULT
      PRINT *, "t_list_prependInt: Invalid link type!"
      CALL sys_halt()
    END SELECT

    ! Integer key
    rlist%Ikey(ipos) = ikey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_prependInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_appendDble(rlist,dkey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends a Double data to the list
!</description>

!<input>
    ! Double key
    REAL(DP), INTENT(IN) :: dkey

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_DOUBLE) THEN
      PRINT *, 't_list_appendDble: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE(rlist%clinkType)
    CASE (LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      END IF
      
    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSE
        rlist%Kprev(ipos)               = rlist%Knext(LTAIL)
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      END IF

    CASE DEFAULT
      PRINT *, "t_list_appendDble: Invalid link type!"
      CALL sys_halt()
    END SELECT

    ! Double key
    rlist%Dkey(ipos)   = dkey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_appendDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_appendSngl(rlist,skey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends a Single data to the list
!</description>

!<input>
    ! Single key
    REAL(SP), INTENT(IN) :: skey

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_SINGLE) THEN
      PRINT *, 't_list_appendSngl: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE(rlist%clinkType)
    CASE (LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      END IF

    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSE
        rlist%Kprev(ipos)               = rlist%Knext(LTAIL)
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      END IF

    CASE DEFAULT
      PRINT *, "t_list_appendSngl: Invalid link type!"
      CALL sys_halt()
    END SELECT

    ! Single key
    rlist%Skey(ipos)   = skey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_appendSngl
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_appendInt(rlist,ikey,ipos,DData,SData,IData)

!<description>
    ! This subroutine appends an Integer data to the list
!</description>

!<input>
    ! Integer key
    INTEGER, INTENT(IN) :: ikey

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the appended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>
    
    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_INT) THEN
      PRINT *, 't_list_appendInt: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE(rlist%clinkType)
    CASE(LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      END IF

    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSE
        rlist%Kprev(ipos)               = rlist%Knext(LTAIL)
        rlist%Knext(rlist%Knext(LTAIL)) = ipos
        rlist%Knext(LTAIL)              = ipos
        rlist%Knext(ipos)               = LNULL
      END IF
      
    CASE DEFAULT
      PRINT *, "t_list_appendInt: Invalid link type!"
      CALL sys_halt()
    END SELECT

    ! Integere key
    rlist%Ikey(ipos)   = ikey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_appendInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_insertDble(rlist,dkey,ipred,ipos,DData,SData,IData)

!<description>
    ! This subroutine inserts a new Double data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Double key
    REAL(DP), INTENT(IN) :: dkey

    ! Position of predecessor
    INTEGER(PREC_LISTIDX), INTENT(IN) :: ipred

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_DOUBLE) THEN
      PRINT *, 't_list_insertDble: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE(rlist%clinkType)
    CASE (LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSEIF (ipred == rlist%Knext(LTAIL)) THEN
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      END IF

    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSEIF (ipred == rlist%Knext(LTAIL)) THEN
        rlist%Kprev(ipos)  = rlist%Knext(LTAIL)
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Kprev(rlist%Knext(ipred)) = ipos
        rlist%Knext(ipos)               = rlist%Knext(ipred)
        rlist%Knext(ipred)              = ipos
        rlist%Kprev(ipos)               = ipred
      END IF

    CASE DEFAULT
      PRINT *, "t_list_insertDble: Invalid link type!"
      CALL sys_halt()
    END SELECT

    ! Double key
    rlist%Dkey(ipos)   = dkey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_insertDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_insertSngl(rlist,skey,ipred,ipos,DData,SData,IData)

!<description>
    ! This subroutine inserts a new Single data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Single key
    REAL(SP), INTENT(IN) :: skey

    ! Position of predecessor
    INTEGER(PREC_LISTIDX), INTENT(IN) :: ipred

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_SINGLE) THEN
      PRINT *, 't_list_insertSngl: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE(rlist%clinkType)
    CASE (LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSEIF (ipred == rlist%Knext(LTAIL)) THEN
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      END IF

    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSEIF (ipred == rlist%Knext(LTAIL)) THEN
        rlist%Kprev(ipos)  = rlist%Knext(LTAIL)
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Kprev(rlist%Knext(ipred)) = ipos
        rlist%Knext(ipos)               = rlist%Knext(ipred)
        rlist%Knext(ipred)              = ipos
        rlist%Kprev(ipos)               = ipred
      END IF

    CASE DEFAULT
      PRINT *, "t_list_insertSngl: Invalid link type!"
      CALL sys_halt()
    END SELECT

    ! Single key
    rlist%Skey(ipos)   = skey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_insertSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_list_insertInt(rlist,ikey,ipred,ipos,DData,IData,SData)

!<description>
    ! This subroutine inserts a new Integer data into the list AFTER
    ! the position ipred
!</description>

!<input>
    ! Integer key
    INTEGER, INTENT(IN) :: ikey

    ! Position of predecessor
    INTEGER(PREC_LISTIDX), INTENT(IN) :: ipred

    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: SData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipos
!</output>
!</subroutine>

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_INT) THEN
      PRINT *, 't_list_insertInt: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    ! Check if list needs to be enlarged
    rlist%na = rlist%na+1
    ipos     = rlist%Knext(LFREE)
    IF (ABS(ipos) > rlist%nna) THEN
      CALL resize(rlist,CEILING(rlist%dfactor*rlist%nna))
    END IF
    
    ! Set next free position
    IF (ipos > 0) THEN
      rlist%Knext(LFREE) = ipos+1
    ELSE
      ipos               = ABS(ipos)
      rlist%Knext(LFREE) = rlist%Knext(ipos)
    END IF
    
    ! Set head, tail and data
    SELECT CASE(rlist%clinkType)
    CASE (LIST_SINGLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSEIF (ipred == rlist%Knext(LTAIL)) THEN
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Knext(ipos)  = rlist%Knext(ipred)
        rlist%Knext(ipred) = ipos
      END IF

    CASE (LIST_DOUBLELINKED)
      IF (rlist%Knext(LHEAD) == LNULL) THEN
        rlist%Knext(LHEAD) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
        rlist%Kprev(ipos)  = LNULL
      ELSEIF (ipred == rlist%Knext(LTAIL)) THEN
        rlist%Kprev(ipos)  = rlist%Knext(LTAIL)
        rlist%Knext(ipred) = ipos
        rlist%Knext(LTAIL) = ipos
        rlist%Knext(ipos)  = LNULL
      ELSE
        rlist%Kprev(rlist%Knext(ipred)) = ipos
        rlist%Knext(ipos)               = rlist%Knext(ipred)
        rlist%Knext(ipred)              = ipos
        rlist%Kprev(ipos)               = ipred
      END IF

    CASE DEFAULT
      PRINT *, "t_list_insertInt: Invalid link type!"
      CALL sys_halt()
    END SELECT

    ! Integer key
    rlist%Ikey(ipos)   = ikey
    
    ! Optional data
    IF ((rlist%isizeInt > 0) .AND. &
        PRESENT(IData)) rlist%IData(:,ipos) = IData
    
    IF ((rlist%isizeDble > 0) .AND. &
        PRESENT(DData)) rlist%DData(:,ipos) = DData
    
    IF ((rlist%isizeSngl > 0) .AND. &
        PRESENT(SData)) rlist%SData(:,ipos) = SData
  END SUBROUTINE t_list_insertInt
  
  ! ***************************************************************************
  
!<function>

  FUNCTION t_list_deleteDble(rlist,dkey) RESULT(f)

!<description>
    ! This function deletes a Double data from the list
!</description>

!<input>
    ! Data
    REAL(DP), INTENT(IN) :: dkey
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_LISTIDX) :: ipred,ipos

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_DOUBLE) THEN
      PRINT *, 't_list_deleteDble: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Search for data
    f=search(rlist,dkey,ipred)
    IF (f == LIST_NOT_FOUND) RETURN
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    IF (rlist%Knext(ipred) == rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  END FUNCTION t_list_deleteDble
  
  ! ***************************************************************************
  
!<function>

  FUNCTION t_list_deleteSngl(rlist,skey) RESULT(f)

!<description>
    ! This function deletes a Single data from the list
!</description>

!<input>
    ! Data
    REAL(SP), INTENT(IN) :: skey
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_LISTIDX) :: ipred,ipos

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_SINGLE) THEN
      PRINT *, 't_list_deleteSngl: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Search for data
    f=search(rlist,skey,ipred)
    IF (f == LIST_NOT_FOUND) RETURN
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    IF (rlist%Knext(ipred) == rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  END FUNCTION t_list_deleteSngl

  ! ***************************************************************************
  
!<function>

  FUNCTION t_list_deleteInt(rlist,ikey) RESULT(f)

!<description>
    ! This function deletes an Integer data from the list
!</description>

!<input>
    ! Data
    INTEGER, INTENT(IN) :: ikey
!</input>

!<inputoutput>
    ! list
    TYPE(t_list), INTENT(INOUT) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_LISTIDX) :: ipred,ipos

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_INT) THEN
      PRINT *, 't_list_deleteInt: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Search for data
    f=search(rlist,ikey,ipred)
    IF (f == LIST_NOT_FOUND) RETURN
    
    ! Delete data
    rlist%na = rlist%na-1
    ipos     = rlist%Knext(ipred)
    IF (rlist%Knext(ipred) == rlist%Knext(LTAIL)) rlist%Knext(LTAIL)=ipred

    ! Update free position
    rlist%Knext(ipred) = rlist%Knext(ipos)
    rlist%Knext(ipos)  = rlist%Knext(LFREE)
    rlist%Knext(LFREE) = -ipos
  END FUNCTION t_list_deleteInt

  ! ***************************************************************************
  
!<function>

  FUNCTION t_list_searchDble(rlist,dkey,ipred) RESULT(f)

!<description>
    ! This function searches for a given Double key in the list
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist

    ! Data
    REAL(DP), INTENT(IN) :: dkey
!</input>

!<output>
    ! Position of the predecessor of the found item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOUT_FOUND / LIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_LISTIDX) :: inext

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_DOUBLE) THEN
      PRINT *, 't_list_searchDble: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)

    ! Check if list is empty
    IF (inext == LNULL) RETURN

    ! What kind of ordering are we
    SELECT CASE(rlist%cordering)
    CASE (LIST_UNORDERED)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%DKey(inext) == dkey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_INCREASING)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%DKey(inext) == dkey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF (rlist%DKey(inext) > dkey) EXIT
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_DECREASING)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%DKey(inext) == dkey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF (rlist%DKey(inext) < dkey) EXIT
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_CSR7)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%DKey(inext) == dkey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF ((ipred.NE.LHEAD) .AND. rlist%DKey(inext) > dkey) EXIT
        
        IF (rlist%Knext(ipred) == rlist%Knext(LTAIL)) THEN
          ipred=rlist%Knext(ipred); EXIT
        END IF
        ipred=rlist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION t_list_searchDble
  
  ! ***************************************************************************
  
!<function>

  FUNCTION t_list_searchSngl(rlist,skey,ipred) RESULT(f)

!<description>
    ! This function searches for a given Single key in the list
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist

    ! Data
    REAL(SP), INTENT(IN) :: skey
!</input>

!<output>
    ! Position of the predecessor of the found item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOT_FOUND / LIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_LISTIDX) :: inext

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_SINGLE) THEN
      PRINT *, 't_list_searchSngl: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)
    
    ! Check if list is empty
    IF (inext == LNULL) RETURN

    ! What kind of ordering are we
    SELECT CASE(rlist%cordering)
    CASE (LIST_UNORDERED)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%SKey(inext) == skey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_INCREASING)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%SKey(inext) == skey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF (rlist%SKey(inext) > skey) EXIT
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_DECREASING)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)  
        IF (rlist%SKey(inext) == skey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF (rlist%SKey(inext) < skey) EXIT
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_CSR7)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%SKey(inext) == skey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF ((ipred.NE.LHEAD) .AND. rlist%SKey(inext) > skey) EXIT
        
        IF (rlist%Knext(ipred) == rlist%Knext(LTAIL)) THEN
          ipred=rlist%Knext(ipred); EXIT
        END IF
        ipred=rlist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION t_list_searchSngl

  ! ***************************************************************************
  
!<function>

  FUNCTION t_list_searchInt(rlist,ikey,ipred) RESULT(f)

!<description>
    ! This function searches for a given Integer key in the list
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist

    ! Data
    INTEGER, INTENT(IN) :: ikey
!</input>

!<output>
    ! Position of the predecessor of the found item
    INTEGER(PREC_LISTIDX), INTENT(OUT) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOT_FOUND / LIST_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_LISTIDX) :: inext

    ! Check if list format is ok
    IF (rlist%clistFormat /= ST_INT) THEN
      PRINT *, 't_list_searchInt: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%Knext(ipred)
    
    ! Check if list is empty
    IF (inext == LNULL) RETURN

    ! What kind of ordering are we
    SELECT CASE(rlist%cordering)
    CASE (LIST_UNORDERED)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%IKey(inext) == ikey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_INCREASING)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%IKey(inext) == ikey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF (rlist%IKey(inext) > ikey) EXIT
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_DECREASING)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%IKey(inext) == ikey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF (rlist%IKey(inext) < ikey) EXIT
        ipred=rlist%Knext(ipred)
      END DO
      
    CASE (LIST_CSR7)
      DO WHILE(ipred.NE.rlist%Knext(LTAIL))
        inext = rlist%Knext(ipred)
        IF (rlist%IKey(inext) == ikey) THEN
          f=LIST_FOUND; EXIT
        END IF
        
        IF ((ipred.NE.LHEAD) .AND. rlist%IKey(inext) > ikey) EXIT
        
        IF (rlist%Knext(ipred) == rlist%Knext(LTAIL)) THEN
          ipred=rlist%Knext(ipred); EXIT
        END IF
        ipred=rlist%Knext(ipred)
      END DO
    END SELECT
  END FUNCTION t_list_searchInt

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE t_list_print(rlist)

!<description>
    ! This subroutine prints the content of the list
!</description>

!<input>
    ! list
    TYPE(t_list), INTENT(IN) :: rlist
!</input>
!</subroutine>

    ! local variable
    INTEGER(PREC_LISTIDX) :: ipos
    
    ipos = rlist%Knext(LHEAD)
    IF (ipos == LNULL) RETURN
    
    SELECT CASE (rlist%clistFormat)
    CASE (ST_DOUBLE)
      DO
        WRITE(*,*) rlist%DKey(ipos)
        IF (ipos == rlist%Knext(LTAIL)) EXIT
        ipos = rlist%Knext(ipos)
      END DO

    CASE (ST_SINGLE)
      DO
        WRITE(*,*) rlist%SKey(ipos)
        IF (ipos == rlist%Knext(LTAIL)) EXIT
        ipos = rlist%Knext(ipos)
      END DO
      
    CASE (ST_INT)
      DO
        WRITE(*,*) rlist%IKey(ipos)
        IF (ipos == rlist%Knext(LTAIL)) EXIT
        ipos = rlist%Knext(ipos)
      END DO
      
    CASE DEFAULT
      PRINT *, 't_list_print: Unsupported data type!'
      CALL sys_halt()
    END SELECT
  END SUBROUTINE t_list_print
END MODULE list
