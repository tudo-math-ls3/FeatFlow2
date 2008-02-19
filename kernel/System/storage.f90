!##############################################################################
!# ****************************************************************************
!# <name> storage </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the FEAT implementation of the global memory
!# management. The usual memory management in Fortran 90 uses pointers
!# for referring to a memory block. As this handling is rather nasty in
!# many circumstances (e.g. it's hard to set up something like an 'array of
!# pointers'), we decided to implement our own memory management - based on
!# ALLOCATE and DEALLOCATE.
!#
!# The new memory management makes use of 'handles'. A handle is an integer
!# valued identifier for a pointer.
!# With 'storage_new', memory can be allocated, a handle is returned.
!# With 'storage_getbase_xxxx', the pointer corresponding to a handle can be
!#  obtained.
!# With 'storage_free', the memory assigned to a handle can be released.
!#
!# The memory management supports 1D and 2D arrays for SINGLE and DOUBLE
!# PRECISION floating point variables as well as I32 integer variables and
!# LOGICAL and CHAR variables.
!#
!# Before any memory is allocated, the memory management must be initialised
!# by a call to 'storage_init'!
!#
!# The following routines can be found here:
!#
!#  1.) storage_init
!#      -> Initialises the storage management
!#
!#  2.) storage_done
!#      -> Cleans up the storage management
!#
!#  3.) storage_info
!#      -> Prints statistics about the heap to the terminal
!#
!#  4.) storage_new  =  storage_new1D / storage_new2D
!#      -> Allocates a new 1D or 2D array
!#
!#  5.) storage_free
!#      -> Releases a handle and the associated memory
!#
!#  6.) storage_getbase_single,
!#      storage_getbase_double,
!#      storage_getbase_int,
!#      storage_getbase_logical,
!#      storage_getbase_char,
!#      -> Determine pointer associated to a handle for singles, doubles,
!#         32-Bit integers, logicals or chars
!#
!#  7.) storage_getbase_single2D,
!#      storage_getbase_double2D,
!#      storage_getbase_int2D,
!#      storage_getbase_logical,
!#      storage_getbase_char,
!#      -> Determine pointer associated to a handle for singles, doubles,
!#         32-Bit integers, logicals or chars, 2D array
!#
!#  8.) storage_copy = storage_copy / storage_copy_explicit / storage_copy_explicit2D
!#      -> Copies the content of one array to another.
!#
!#  9.) storage_clear
!#      -> Clears an array by overwriting the entries with 0 (or .FALSE. for LOGICALs).
!#
!# 10.) storage_getsize = storage_getsize1d / storage_getsize2d
!#      -> Get the length of an array on the heap.
!#
!# 11.) storage_getdatatype
!#      -> Get the datatype of an array on the heap.
!#
!# 12.) storage_getdimension
!#      -> Get the dimension of an array on the heap.
!#
!# 13.) storage_realloc
!#      -> Reallocate a 1D or 2D array (only 2nd dimension)
!#
!# 14.) storage_initialiseBlock
!#      -> Initialise a storage block with zero (like storage_clear)
!#         or increasing number.
!#
!# 15.) storage_isEqual
!#      -> Checks if the content of two different handles is equal
!# </purpose>
!##############################################################################

MODULE storage

  USE fsystem
  USE linearalgebra
  USE genoutput

  IMPLICIT NONE

!<constants>

  !<constantblock description="Storage block type identifiers">

  ! defines an non-allocated storage block handle
  INTEGER, PARAMETER :: ST_NOHANDLE = 0

  ! storage block contains single floats
  INTEGER, PARAMETER :: ST_SINGLE = 1

  ! storage block contains double floats
  INTEGER, PARAMETER :: ST_DOUBLE = 2

  ! storage block contains ints
  INTEGER, PARAMETER :: ST_INT = 3

  ! storage block contains logicals
  INTEGER, PARAMETER :: ST_LOGICAL = 4

  ! storage block contains characters
  INTEGER, PARAMETER :: ST_CHAR = 5

  !</constantblock>

  !<constantblock description="Constants for initialisation of memory on allocation">

  ! init new storage block with zeros (or .FALSE. for logicals)
  INTEGER, PARAMETER :: ST_NEWBLOCK_ZERO = 0

  ! no init new storage block
  INTEGER, PARAMETER :: ST_NEWBLOCK_NOINIT = 1

  ! init new storage block with 1,2,3,...,n
  ! Note: This initialization constant must NOT be used for initialisation of logicals
  !       or characters!!!
  INTEGER, PARAMETER :: ST_NEWBLOCK_ORDERED = 2

  !</constantblock>

  !<constantblock description="Constants for calculating memory">

  ! How many bytes has an integer?
  INTEGER :: ST_INT2BYTES = I32

  ! How many bytes has a standard real?
  INTEGER :: ST_SINGLE2BYTES = SP

  ! How many bytes has a double precision real?
  INTEGER :: ST_DOUBLE2BYTES = DP

  ! How many bytes has a logical?
  INTEGER :: ST_LOGICAL2BYTES = I32

  ! How many bytes has a character?
  ! Note: We are not 100% sure, but this may differ on other architectures... O_o
  INTEGER :: ST_CHAR2BYTES = 1

  !</constantblock>

!</constants>

!<types>

  !<typeblock>

  ! Type block for describing a handle. This collects the number of the
  ! handle, the storage amount associated to it, the pointer to the memory
  ! location etc.

  TYPE t_storageNode
    PRIVATE

    ! Type of data associated to the handle (ST_NOHANDLE, ST_SINGLE,
    ! ST_DOUBLE, ST_INT)
    INTEGER :: idataType = ST_NOHANDLE

    ! Dimension associated to the handle (0=not assigned, 1=1D, 2=2D array)
    INTEGER :: idimension = 0

    ! The name of the array that is associated to that handle
    CHARACTER(LEN=SYS_NAMELEN) :: sname = ''

    ! Amount of memory (in bytes) associated to this block.
    ! We store that as a double to allow storing numbers > 2GB !
    REAL(DP) :: dmemBytes = 0.0_DP

    ! Pointer to 1D real array or NULL() if not assigned
    REAL(SP), DIMENSION(:), POINTER       :: p_Fsingle1D    => NULL()

    ! Pointer to 1D double precision array or NULL() if not assigned
    REAL(DP), DIMENSION(:), POINTER       :: p_Ddouble1D    => NULL()

    ! Pointer to 1D integer array or NULL() if not assigned
    INTEGER(I32), DIMENSION(:), POINTER   :: p_Iinteger1D   => NULL()

    ! Pointer to 1D logical array or NULL() if not assigned
    LOGICAL, DIMENSION(:), POINTER        :: p_Blogical1D   => NULL()

    ! Pointer to 1D character array or NULL() if not assigned
    CHARACTER, DIMENSION(:), POINTER      :: p_Schar1D      => NULL()

    ! Pointer to 2D real array or NULL() if not assigned
    REAL(SP), DIMENSION(:,:), POINTER     :: p_Fsingle2D    => NULL()

    ! Pointer to 2D double precision array or NULL() if not assigned
    REAL(DP), DIMENSION(:,:), POINTER     :: p_Ddouble2D    => NULL()

    ! Pointer to 2D integer array or NULL() if not assigned
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_Iinteger2D   => NULL()

    ! Pointer to 2D logical array or NULL() if not assigned
    LOGICAL, DIMENSION(:,:), POINTER        :: p_Blogical2D => NULL()

    ! Pointer to 2D character array or NULL() if not assigned
    CHARACTER, DIMENSION(:,:), POINTER      :: p_Schar2D    => NULL()

  END TYPE

  !</typeblock>

  !<typeblock>

  ! This block represents a heap that maintains singole, double precision
  ! and integer data. It contains a list of t_storageNode elements for all
  ! the handles.
  ! There's one global object of this type for the global storage management,
  ! but if necessary, an algorithm can create such a block locally, too,
  ! to prevent conflicts with the global memory.

  TYPE t_storageBlock

    PRIVATE

    ! An array of t_storageNode objects corresponding to the handles.
    ! Can be dynamically extended if there are not enough handles available.
    TYPE(t_storageNode), DIMENSION(:), POINTER :: p_Rdescriptors => NULL()

    ! A list of all 'free' handles. This is a 'ring' queue. If all
    ! handles are in use, p_Rdescriptors and p_IfreeHandles are dynamically
    ! extended.
    INTEGER, DIMENSION(:), POINTER :: p_IfreeHandles => NULL()

    ! Index in p_IfreeHandles to the next free handle
    INTEGER :: p_inextFreeHandle = 0

    ! Index in p_IfreeHandles to the last free handle
    INTEGER :: p_ilastFreeHandle = 0

    ! Number of handles in use
    INTEGER :: ihandlesInUse = 0

    ! Total number of handles maintained by this block; = size(p_Rdescriptors).
    INTEGER :: nhandlesTotal = 0

    ! Number of handles to add if there are not enough free handles.
    INTEGER :: ihandlesDelta = 0

    ! Total amount of memory (in bytes) that is in use. We maintain it
    ! as a double as this allows to save values > 2GB!
    REAL(DP) :: dtotalMem = 0.0_DP

    ! Maximum number of handles that were in use ofer the whole lifetime
    ! of this structure.
    INTEGER :: nhandlesInUseMax = 0

    ! Maximum amount of memory that was in use over the whole lifetime
    ! of this structure.
    REAL(DP) :: dtotalMemMax = 0.0_DP

  END TYPE

  !</typeblock>

!</types>

!<globals>

  ! Global memory management structure
  TYPE(t_storageBlock), SAVE, TARGET :: rbase

!</globals>

  INTERFACE storage_new
    MODULE PROCEDURE storage_new1D
    MODULE PROCEDURE storage_new1DFixed
    MODULE PROCEDURE storage_new2D
    MODULE PROCEDURE storage_new2DFixed
  END INTERFACE

  INTERFACE storage_realloc
    MODULE PROCEDURE storage_realloc
    MODULE PROCEDURE storage_reallocFixed
  END INTERFACE
  
  INTERFACE storage_getbase_int
    MODULE PROCEDURE storage_getbase_int
    MODULE PROCEDURE storage_getbase_intUBnd
    MODULE PROCEDURE storage_getbase_intLUBnd
  END INTERFACE

  INTERFACE storage_getbase_single
    MODULE PROCEDURE storage_getbase_single
    MODULE PROCEDURE storage_getbase_singleUBnd
    MODULE PROCEDURE storage_getbase_singleLUBnd
  END INTERFACE

  INTERFACE storage_getbase_double
    MODULE PROCEDURE storage_getbase_double
    MODULE PROCEDURE storage_getbase_doubleUBnd
    MODULE PROCEDURE storage_getbase_doubleLUBnd
  END INTERFACE

  INTERFACE storage_getbase_logical
    MODULE PROCEDURE storage_getbase_logical
    MODULE PROCEDURE storage_getbase_logicalUBnd
    MODULE PROCEDURE storage_getbase_logicalLUBnd
  END INTERFACE

  INTERFACE storage_getbase_char
    MODULE PROCEDURE storage_getbase_char
    MODULE PROCEDURE storage_getbase_charUBnd
    MODULE PROCEDURE storage_getbase_charLUBnd
  END INTERFACE

  INTERFACE storage_getbase_int2D
    MODULE PROCEDURE storage_getbase_int2D
    MODULE PROCEDURE storage_getbase_int2DUBnd
    MODULE PROCEDURE storage_getbase_int2DLUBnd
  END INTERFACE

  INTERFACE storage_getbase_single2D
    MODULE PROCEDURE storage_getbase_single2D
    MODULE PROCEDURE storage_getbase_single2DUBnd
    MODULE PROCEDURE storage_getbase_single2DLUBnd
  END INTERFACE

  INTERFACE storage_getbase_double2D
    MODULE PROCEDURE storage_getbase_double2D
    MODULE PROCEDURE storage_getbase_double2DUBnd
    MODULE PROCEDURE storage_getbase_double2DLUBnd
  END INTERFACE

  INTERFACE storage_getbase_logical2D
    MODULE PROCEDURE storage_getbase_logical2D
    MODULE PROCEDURE storage_getbase_logical2DUBnd
    MODULE PROCEDURE storage_getbase_logical2DLUBnd
  END INTERFACE

  INTERFACE storage_getbase_char2D
    MODULE PROCEDURE storage_getbase_char2D
    MODULE PROCEDURE storage_getbase_char2DUBnd
    MODULE PROCEDURE storage_getbase_char2DLUBnd
  END INTERFACE

  INTERFACE storage_getsize
    MODULE PROCEDURE storage_getsize1D
    MODULE PROCEDURE storage_getsize2D
  END INTERFACE

  INTERFACE storage_copy
    MODULE PROCEDURE storage_copy
    MODULE PROCEDURE storage_copy_explicit
    MODULE PROCEDURE storage_copy_explicit2D
  END INTERFACE

CONTAINS

!************************************************************************

!<subroutine>

  SUBROUTINE storage_init(ihandleCount, ihandlesDelta, rheap)

!<description>

  ! This routine initializes the storage management.
  ! ihandleCount is the initial number of handles maintained by the
  ! storage routines. If there are not enough free handles, the number
  ! of handles are increased by ihandlesDelta (which is initially set
  ! to 1/2*ihandleCount if not given).
  ! rheap allows to specify a 'local' heap structure to initialise.
  ! If not given, the global memory management is initialised.

!</description>

!<input>

  ! Initial number of handles maintained by the storage routines.
  INTEGER, INTENT(IN) :: ihandleCount

  ! OPTIONAL: Number of handles to increase the memory block by, if there are
  ! not enough handles available. Standard setting is 1/2*ihandleCount.
  INTEGER, INTENT(IN), OPTIONAL :: ihandlesDelta

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is initialised.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!</subroutine>

  ! local variables

  ! the real 'handle-delta'
  INTEGER :: ihandles, ihDelta

  ! Pointer to the heap to initialise
  TYPE(t_storageBlock), POINTER :: p_rheap

  INTEGER :: i

  ! Initialise ihDelta and p_rheap and work with these - as the other
  ! parameters are optional.
  ! We work at least with 1 handles and ihDelta = 1.

  ihandles = MAX(1,ihandlecount)

  ihDelta = 1
  IF(PRESENT(ihandlesDelta)) ihDelta = ihandlesDelta
  ihDelta = MAX(1,ihDelta)

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  ! Initialise the memory management block

  p_rheap%nhandlesTotal = ihandles
  p_rheap%ihandlesDelta = ihDelta
  p_rheap%p_inextFreeHandle = 1
  p_rheap%p_ilastFreeHandle = ihandles
  p_rheap%ihandlesInUse = 0
  p_rheap%nhandlesInUseMax = 0
  ALLOCATE(p_rheap%p_Rdescriptors(ihandles))
  ALLOCATE(p_rheap%p_IfreeHandles(ihandles))

  ! All handles free
  DO i=1,ihandles
    p_rheap%p_IfreeHandles(i) = i
  END DO

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_done(rheap)

!<description>
  ! This routine cleans up the storage management. All data on the
  ! heap is released from memory.
!</description>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is cleaned up.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap to initialise
  TYPE(t_storageBlock), POINTER :: p_rheap

  INTEGER :: i,ihandle

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  ! Delete all data from the heap
  DO i = 1,SIZE(p_rheap%p_Rdescriptors)
    ! Don't pass i as handle as storage_free will set the handle
    ! passed to it to 0!
    ihandle = i
    IF (p_rheap%p_Rdescriptors(i)%idataType .NE. ST_NOHANDLE) &
      CALL storage_free(ihandle,rheap)
  END DO

  ! Clean up the memory management block
  p_rheap%nhandlesTotal = 0
  p_rheap%ihandlesDelta = 0
  p_rheap%p_inextFreeHandle = 0
  p_rheap%p_ilastFreeHandle = 0
  p_rheap%ihandlesInUse = 0

  ! Release the descriptors
  DEALLOCATE(p_rheap%p_IfreeHandles)
  DEALLOCATE(p_rheap%p_Rdescriptors)

  END SUBROUTINE

!************************************************************************

!<function>

  INTEGER FUNCTION storage_newhandle (rheap) RESULT(ihandle)

!<description>
  ! This routine creates a new handle in the heap structure rheap and
  ! returns the handle number in ihandle. If there is no handle
  ! available, the heap structure is changed to take more handles.
!</description>

!<result>
  ! The new handle number.
!</result>

!<inputoutput>

  ! The heap structure where to create a new handle
  TYPE(t_storageBlock), INTENT(INOUT) :: rheap

!</inputoutput>

!</function>

  ! local variables
  TYPE(t_storageNode), DIMENSION(:), POINTER :: p_Rdescriptors => NULL()
  INTEGER, DIMENSION(:), POINTER :: p_IfreeHandles => NULL()
  INTEGER :: i

  IF (rheap%nhandlesTotal .LE. 0) THEN
    PRINT *,'Error: Heap not initialised!'
    CALL sys_halt()
  END IF

  ! Handles available?

  IF (rheap%ihandlesInUse .GE. rheap%nhandlesTotal) THEN

    ! All handles are in use. We have to modify our ring to accept more
    ! handles.
    !
    ! At first, reallocate the descriptor-array and the queue-array with
    ! the new size.
    ALLOCATE (p_Rdescriptors (rheap%nhandlesTotal + rheap%ihandlesDelta) )
    ALLOCATE (p_IfreeHandles (rheap%nhandlesTotal + rheap%ihandlesDelta) )

    ! Copy the content, release the old arrays and replace them by the new
    ! ones.
    p_Rdescriptors(1:rheap%nhandlesTotal) = rheap%p_Rdescriptors(1:rheap%nhandlesTotal)
    p_IfreeHandles(1:rheap%nhandlesTotal) = rheap%p_IfreeHandles(1:rheap%nhandlesTotal)

    DEALLOCATE(rheap%p_Rdescriptors)
    DEALLOCATE(rheap%p_IfreeHandles)
    rheap%p_Rdescriptors => p_Rdescriptors
    rheap%p_IfreeHandles => p_IfreeHandles

    ! Add the new handles to the list of 'free' handles.
    DO i=rheap%nhandlesTotal+1, rheap%nhandlesTotal + rheap%ihandlesDelta
      p_IfreeHandles (i) = i
    END DO

    ! The first new 'free' handle is not at position...
    rheap%p_inextFreeHandle = rheap%nhandlesTotal+1

    ! And the last 'free' handle is at the end of the new list.
    rheap%p_ilastFreeHandle = rheap%nhandlesTotal + rheap%ihandlesDelta

    ! Modify the heap structure - we have more handles now.
    rheap%nhandlesTotal = rheap%nhandlesTotal + rheap%ihandlesDelta

  END IF

  ! Get the new handle...
  ihandle = rheap%p_IfreeHandles (rheap%p_inextFreeHandle)

  ! and modify our queue pointers that we use a new one.
  rheap%p_inextFreeHandle = MOD(rheap%p_inextFreeHandle,rheap%nhandlesTotal)+1

  rheap%ihandlesInUse = rheap%ihandlesInUse + 1

  rheap%nhandlesInUseMax = MAX(rheap%nhandlesInUseMax,rheap%ihandlesInUse)

  END FUNCTION

!************************************************************************

!<subroutine>

  SUBROUTINE storage_releasehandle (ihandle,rheap)

!<description>
  ! This routine releases a handle from the heap structure rheap.
  ! Memory is not deallocated, simply the structures are cleaned up.
!</description>

!<input>
  ! The handle to release
  INTEGER, INTENT(INOUT) :: ihandle
!</input>

!<inputoutput>
  ! The heap structure where to release the handle from.
  TYPE(t_storageBlock), INTENT(INOUT) :: rheap
!</inputoutput>

!</subroutine>

  TYPE(t_storageNode), POINTER :: p_rnode

  ! Where is the descriptor of the handle?
  p_rnode => rheap%p_Rdescriptors(ihandle)

  ! Subtract the memory amount from the statistics
  rheap%dtotalMem = rheap%dtotalMem - p_rnode%dmemBytes

  ! Clear the descriptor structure
  p_rnode%idataType = ST_NOHANDLE
  p_rnode%idimension = 0
  p_rnode%dmemBytes = 0.0_DP
  NULLIFY(p_rnode%p_Fsingle1D)
  NULLIFY(p_rnode%p_Ddouble1D)
  NULLIFY(p_rnode%p_Iinteger1D)
  NULLIFY(p_rnode%p_Blogical1D)
  NULLIFY(p_rnode%p_Schar1D)
  NULLIFY(p_rnode%p_Fsingle2D)
  NULLIFY(p_rnode%p_Ddouble2D)
  NULLIFY(p_rnode%p_Iinteger2D)
  NULLIFY(p_rnode%p_Blogical2D)
  NULLIFY(p_rnode%p_Schar2D)

  ! Handle ihandle is available now - put it to the list of available handles.
  rheap%p_ilastFreeHandle = MOD(rheap%p_ilastFreeHandle,rheap%nhandlesTotal) + 1
  rheap%p_IfreeHandles (rheap%p_ilastFreeHandle) = ihandle

  rheap%ihandlesInUse = rheap%ihandlesInUse - 1

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_initialiseNode (rstorageNode,cinitNewBlock,istartIndex,istopIndex)

!<description>
  ! Internal subroutine: Initialise the memory identified by storage
  ! node rstorageNode according to the constant cinitNewBlock.
!</description>

!<input>
  ! The storage node whose associated storage should be initialised.
  TYPE(t_storageNode), INTENT(INOUT) :: rstorageNode

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO, ST_NEWBLOCK_NOINIT,
  ! ST_NEWBLOCK_ORDERED). Specifies how to initialise the data block associated
  ! to ihandle.
  INTEGER, INTENT(IN) :: cinitNewBlock

  ! Start index from where to initialise; should usually be =1.
  ! For multidimensional arrays, this specifies the start index of the
  ! last dimension.
  INTEGER(I32), INTENT(IN) :: istartIndex

  ! OPTIONAL: Stop index up to which to initialise; should usually be =SIZE(*).
  ! For multidimensional arrays, this specifies the stop index of the
  ! last dimension.
  INTEGER(I32), INTENT(IN), OPTIONAL :: istopIndex
!</input>

!</subroutine>

    ! variable for ordering 1,2,3,...,N
    INTEGER :: iorder

    SELECT CASE (rstorageNode%idimension)
    CASE (1)

      SELECT CASE (cinitNewBlock)
      CASE (ST_NEWBLOCK_ZERO)
        ! Clear the vector if necessary
        IF (PRESENT(istopIndex)) THEN
          SELECT CASE (rstorageNode%idataType)
          CASE (ST_SINGLE)
            rstorageNode%p_Fsingle1D(istartIndex:istopIndex) = 0.0_SP
          CASE (ST_DOUBLE)
            rstorageNode%p_Ddouble1D(istartIndex:istopIndex) = 0.0_DP
          CASE (ST_INT)
            rstorageNode%p_Iinteger1D(istartIndex:istopIndex) = 0_I32
          CASE (ST_LOGICAL)
            rstorageNode%p_Blogical1D(istartIndex:istopIndex) = .FALSE.
          CASE (ST_CHAR)
            rstorageNode%p_Schar1D(istartIndex:istopIndex) = ACHAR(0)
          END SELECT
        ELSE
          SELECT CASE (rstorageNode%idataType)
          CASE (ST_SINGLE)
            rstorageNode%p_Fsingle1D(istartIndex:) = 0.0_SP
          CASE (ST_DOUBLE)
            rstorageNode%p_Ddouble1D(istartIndex:) = 0.0_DP
          CASE (ST_INT)
            rstorageNode%p_Iinteger1D(istartIndex:) = 0_I32
          CASE (ST_LOGICAL)
            rstorageNode%p_Blogical1D(istartIndex:) = .FALSE.
          CASE (ST_CHAR)
            rstorageNode%p_Schar1D(istartIndex:) = ACHAR(0)
          END SELECT
        END IF

      CASE (ST_NEWBLOCK_ORDERED)
        ! Impose ordering 1,2,3,...,N if necessary
        IF (PRESENT(istopIndex)) THEN
          SELECT CASE (rstorageNode%idataType)
          CASE (ST_SINGLE)
            DO iorder=istartIndex,istopIndex
              rstorageNode%p_Fsingle1D(iorder) = REAL(iorder,SP)
            END DO
          CASE (ST_DOUBLE)
            DO iorder=istartIndex,istopIndex
              rstorageNode%p_Ddouble1D(iorder) = REAL(iorder,DP)
            END DO
          CASE (ST_INT)
            DO iorder=istartIndex,istopIndex
              rstorageNode%p_Iinteger1D(iorder) = INT(iorder,I32)
            END DO
          CASE (ST_LOGICAL)
            PRINT *, "Error: Logical array can not be initialised with ST_NEWBLOCK_ORDERED !"
            CALL sys_halt()
          CASE (ST_CHAR)
            PRINT *, "Error: Character array can not be initialised with ST_NEWBLOCK_ORDERED !"
            CALL sys_halt()
          END SELECT
        ELSE
          SELECT CASE (rstorageNode%idataType)
          CASE (ST_SINGLE)
            DO iorder=istartIndex,UBOUND(rstorageNode%p_Fsingle1D,1)
              rstorageNode%p_Fsingle1D(iorder) = REAL(iorder,SP)
            END DO
          CASE (ST_DOUBLE)
            DO iorder=istartIndex,UBOUND(rstorageNode%p_Ddouble1D,1)
              rstorageNode%p_Ddouble1D(iorder) = REAL(iorder,DP)
            END DO
          CASE (ST_INT)
            DO iorder=istartIndex,UBOUND(rstorageNode%p_Iinteger1D,1)
              rstorageNode%p_Iinteger1D(iorder) = INT(iorder,I32)
            END DO
          CASE (ST_LOGICAL)
            PRINT *, "Error: Logical array can not be initialised with ST_NEWBLOCK_ORDERED !"
            CALL sys_halt()
          CASE (ST_CHAR)
            PRINT *, "Error: Character array can not be initialised with ST_NEWBLOCK_ORDERED !"
            CALL sys_halt()
          END SELECT
        END IF

      END SELECT

    CASE (2)

      SELECT CASE (cinitNewBlock)
      CASE (ST_NEWBLOCK_ZERO)
        ! Clear the vector
        IF (PRESENT(istopIndex)) THEN
          SELECT CASE (rstorageNode%idataType)
          CASE (ST_SINGLE)
            rstorageNode%p_Fsingle2D(:,istartIndex:istopIndex) = 0.0_SP
          CASE (ST_DOUBLE)
            rstorageNode%p_Ddouble2D(:,istartIndex:istopIndex) = 0.0_DP
          CASE (ST_INT)
            rstorageNode%p_Iinteger2D(:,istartIndex:istopIndex) = 0_I32
          CASE (ST_LOGICAL)
            rstorageNode%p_Blogical2D(:,istartIndex:istopIndex) = .FALSE.
          CASE (ST_CHAR)
            rstorageNode%p_Schar2D(:,istartIndex:istopIndex) = ACHAR(0)
          END SELECT
        ELSE
          SELECT CASE (rstorageNode%idataType)
          CASE (ST_SINGLE)
            rstorageNode%p_Fsingle2D(:,istartIndex:) = 0.0_SP
          CASE (ST_DOUBLE)
            rstorageNode%p_Ddouble2D(:,istartIndex:) = 0.0_DP
          CASE (ST_INT)
            rstorageNode%p_Iinteger2D(:,istartIndex:) = 0_I32
          CASE (ST_LOGICAL)
            rstorageNode%p_Blogical2D(:,istartIndex:) = .FALSE.
          CASE (ST_CHAR)
            rstorageNode%p_Schar2D(:,istartIndex:) = ACHAR(0)
          END SELECT
        END IF

      CASE (ST_NEWBLOCK_ORDERED)
        PRINT *, 'Error: ordering not available for multidimensional array'
        CALL sys_halt()

      END SELECT

    CASE DEFAULT
      PRINT *,'storage_initialiseNode: Unsupported dimension'
      CALL sys_halt()

    END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_initialiseBlock (ihandle, cinitNewBlock,&
      rheap, istartIndex)

!<description>
  ! This routine initialises the memory associated to ihandle according
  ! to the constant cinitNewBlock.
!</description>

!<input>
  ! Handle of the memory block to initialise
  INTEGER, INTENT(IN) :: ihandle

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO, ST_NEWBLOCK_NOINIT,
  ! ST_NEWBLOCK_ORDERED). Specifies how to initialise the data block associated
  ! to ihandle.
  INTEGER, INTENT(IN) :: cinitNewBlock

  ! OPTIONAL: Start index of Block
  INTEGER(I32), INTENT(IN), OPTIONAL :: istartIndex
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!</subroutine>

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Initialise
  IF (PRESENT(istartIndex)) THEN
    CALL storage_initialiseNode (p_rnode,cinitNewBlock,istartIndex)
  ELSE
    CALL storage_initialiseNode (p_rnode,cinitNewBlock,1_I32)
  END IF

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_new1D (scall, sname, isize, ctype, ihandle, &
                            cinitNewBlock, rheap)

!<description>
  !This routine reserves a 1D memory block of desired size and type.
!</description>

!<input>

  !name of the calling routine
  CHARACTER(LEN=*), INTENT(IN) :: scall

  !clear name of data field
  CHARACTER(LEN=*), INTENT(IN) :: sname

  !requested storage size
  INTEGER(I32), INTENT(IN) :: isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  INTEGER, INTENT(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT,ST_NEWBLOCK_ORDERED)
  INTEGER, INTENT(IN) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  INTEGER, INTENT(OUT) :: ihandle

!</output>

!</subroutine>

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode
  CHARACTER(LEN=SYS_NAMELEN) :: snameBackup

  IF (isize .EQ. 0) THEN
    PRINT *,'storage_new1D Warning: isize=0'
    ihandle = ST_NOHANDLE
    RETURN
  END IF

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF
  
  ! Back up the array name. This is important for a very crual situation:
  ! The storage_newhandle below may reallocate the p_Rdescriptors array.
  ! If sname points to one of the old arrays, the pointer gets invalid
  ! and the name cannot be accessed anymore. So make a backup of that 
  ! before creating a new handle!
  snameBackup = sname

  ! Get a new handle.
  ihandle = storage_newhandle (p_rheap)

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Initialise the content

  p_rnode%idataType = ctype
  p_rnode%idimension = 1
  p_rnode%sname = snameBackup

  ! Allocate memory according to isize:

  SELECT CASE (ctype)
  CASE (ST_SINGLE)
    ALLOCATE(p_rnode%p_Fsingle1D(isize))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_SINGLE2BYTES)
  CASE (ST_DOUBLE)
    ALLOCATE(p_rnode%p_Ddouble1D(isize))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_DOUBLE2BYTES)
  CASE (ST_INT)
    ALLOCATE(p_rnode%p_Iinteger1D(isize))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_INT2BYTES)
  CASE (ST_LOGICAL)
    ALLOCATE(p_rnode%p_Blogical1D(isize))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_LOGICAL2BYTES)
  CASE (ST_CHAR)
    ALLOCATE(p_rnode%p_Schar1D(isize))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_CHAR2BYTES)
  CASE DEFAULT
    PRINT *,'Error: unknown mem type'
    CALL sys_halt()
  END SELECT

  p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
  IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
    p_rheap%dtotalMemMax = p_rheap%dtotalMem

  ! Initialise the memory block
  CALL storage_initialiseBlock (ihandle, cinitNewBlock, rheap)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_new1Dfixed (scall, sname, ilbound, iubound,&
                            ctype, ihandle, cinitNewBlock, rheap)

!<description>
  !This routine reserves a 1D memory block of desired bounds and type.
!</description>

!<input>

  !name of the calling routine
  CHARACTER(LEN=*), INTENT(IN) :: scall

  !clear name of data field
  CHARACTER(LEN=*), INTENT(IN) :: sname

  !requested lower bound
  INTEGER(I32), INTENT(IN) :: ilbound

  !requested upper bound
  INTEGER(I32), INTENT(IN) :: iubound

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  INTEGER, INTENT(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT,ST_NEWBLOCK_ORDERED)
  INTEGER, INTENT(IN) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  INTEGER, INTENT(OUT) :: ihandle

!</output>

!</subroutine>

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode
  INTEGER(I32) :: isize
  CHARACTER(LEN=SYS_NAMELEN) :: snameBackup

  isize=iubound-ilbound+1
  IF (isize .EQ. 0) THEN
    PRINT *,'storage_new1Dfixed Warning: isize=0'
    ihandle = ST_NOHANDLE
    RETURN
  END IF

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  ! Back up the array name. This is important for a very crual situation:
  ! The storage_newhandle below may reallocate the p_Rdescriptors array.
  ! If sname points to one of the old arrays, the pointer gets invalid
  ! and the name cannot be accessed anymore. So make a backup of that 
  ! before creating a new handle!
  snameBackup = sname

  ! Get a new handle
  ihandle = storage_newhandle (p_rheap)

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Initialise the content

  p_rnode%idataType = ctype
  p_rnode%idimension = 1
  p_rnode%sname = snameBackup

  ! Allocate memory according to isize:

  SELECT CASE (ctype)
  CASE (ST_SINGLE)
    ALLOCATE(p_rnode%p_Fsingle1D(ilbound:iubound))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_SINGLE2BYTES)
  CASE (ST_DOUBLE)
    ALLOCATE(p_rnode%p_Ddouble1D(ilbound:iubound))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_DOUBLE2BYTES)
  CASE (ST_INT)
    ALLOCATE(p_rnode%p_Iinteger1D(ilbound:iubound))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_INT2BYTES)
  CASE (ST_LOGICAL)
    ALLOCATE(p_rnode%p_Blogical1D(ilbound:iubound))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_LOGICAL2BYTES)
  CASE (ST_CHAR)
    ALLOCATE(p_rnode%p_Schar1D(ilbound:iubound))
    p_rnode%dmemBytes = REAL(isize,DP)*REAL(ST_CHAR2BYTES)
  CASE DEFAULT
    PRINT *,'Error: unknown mem type'
    CALL sys_halt()
  END SELECT

  p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
  IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
    p_rheap%dtotalMemMax = p_rheap%dtotalMem

  ! Initialise the memory block
  CALL storage_initialiseBlock (ihandle, cinitNewBlock, rheap, ilbound)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_new2D (scall, sname, Isize, ctype, ihandle, &
                            cinitNewBlock, rheap)

!<description>
  !This routine reserves a 2D memory block of desired size and type.
!</description>

!<input>

  !name of the calling routine
  CHARACTER(LEN=*), INTENT(IN) :: scall

  !clear name of data field
  CHARACTER(LEN=*), INTENT(IN) :: sname

  !requested storage size for 1st and 2nd dimension
  INTEGER(I32), DIMENSION(2), INTENT(IN) :: Isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  INTEGER, INTENT(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  INTEGER, INTENT(IN) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  INTEGER, INTENT(OUT) :: ihandle

!</output>

!</subroutine>

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode
  CHARACTER(LEN=SYS_NAMELEN) :: snameBackup

  IF ((Isize(1) .EQ. 0) .OR. (Isize(2) .EQ. 0)) THEN
    PRINT *,'storage_new2D Warning: Isize=0'
    ihandle = 0
    RETURN
  END IF

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  ! Back up the array name. This is important for a very crual situation:
  ! The storage_newhandle below may reallocate the p_Rdescriptors array.
  ! If sname points to one of the old arrays, the pointer gets invalid
  ! and the name cannot be accessed anymore. So make a backup of that 
  ! before creating a new handle!
  snameBackup = sname

  ! Get a new handle
  ihandle = storage_newhandle (p_rheap)

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Initialise the content

  p_rnode%idataType = ctype
  p_rnode%idimension = 2
  p_rnode%sname = snameBackup

  ! Allocate memory according to Isize:

  SELECT CASE (ctype)
  CASE (ST_SINGLE)
    ALLOCATE(p_rnode%p_Fsingle2D(Isize(1),Isize(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_SINGLE2BYTES)
  CASE (ST_DOUBLE)
    ALLOCATE(p_rnode%p_Ddouble2D(Isize(1),Isize(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_DOUBLE2BYTES)
  CASE (ST_INT)
    ALLOCATE(p_rnode%p_Iinteger2D(Isize(1),Isize(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_INT2BYTES)
  CASE (ST_LOGICAL)
    ALLOCATE(p_rnode%p_Blogical2D(Isize(1),Isize(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_LOGICAL2BYTES)
  CASE (ST_CHAR)
    ALLOCATE(p_rnode%p_Schar2D(Isize(1),Isize(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_CHAR2BYTES)
  CASE DEFAULT
    PRINT *,'Error: unknown mem type'
    CALL sys_halt()
  END SELECT

  p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
  IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
    p_rheap%dtotalMemMax = p_rheap%dtotalMem

  ! Initialise the storage block
  CALL storage_initialiseBlock (ihandle, cinitNewBlock, rheap)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_new2Dfixed (scall, sname, Ilbound, Iubound, ctype,&
                            ihandle, cinitNewBlock, rheap)

!<description>
  !This routine reserves a 2D memory block of desired bounds and type.
!</description>

!<input>

  !name of the calling routine
  CHARACTER(LEN=*), INTENT(IN) :: scall

  !clear name of data field
  CHARACTER(LEN=*), INTENT(IN) :: sname

  !requested lower bounds for 1st and 2nd dimension
  INTEGER(I32), DIMENSION(2), INTENT(IN) :: Ilbound

  !requested upper bounds for 1st and 2nd dimension
  INTEGER(I32), DIMENSION(2), INTENT(IN) :: Iubound

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  INTEGER, INTENT(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  INTEGER, INTENT(IN) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  INTEGER, INTENT(OUT) :: ihandle

!</output>

!</subroutine>

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode
  INTEGER(I32), DIMENSION(2) :: Isize
  CHARACTER(LEN=SYS_NAMELEN) :: snameBackup

  Isize=Iubound-Ilbound+1
  IF ((Isize(1) .EQ. 0) .OR. (Isize(2) .EQ. 0)) THEN
    PRINT *,'storage_new2D Warning: Isize=0'
    ihandle = 0
    RETURN
  END IF

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  ! Back up the array name. This is important for a very crual situation:
  ! The storage_newhandle below may reallocate the p_Rdescriptors array.
  ! If sname points to one of the old arrays, the pointer gets invalid
  ! and the name cannot be accessed anymore. So make a backup of that 
  ! before creating a new handle!
  snameBackup = sname

  ! Get a new handle
  ihandle = storage_newhandle (p_rheap)

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Initialise the content

  p_rnode%idataType = ctype
  p_rnode%idimension = 2
  p_rnode%sname = snameBackup

  ! Allocate memory according to Isize:

  SELECT CASE (ctype)
  CASE (ST_SINGLE)
    ALLOCATE(p_rnode%p_Fsingle2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_SINGLE2BYTES)
  CASE (ST_DOUBLE)
    ALLOCATE(p_rnode%p_Ddouble2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_DOUBLE2BYTES)
  CASE (ST_INT)
    ALLOCATE(p_rnode%p_Iinteger2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_INT2BYTES)
  CASE (ST_LOGICAL)
    ALLOCATE(p_rnode%p_Blogical2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_LOGICAL2BYTES)
  CASE (ST_CHAR)
    ALLOCATE(p_rnode%p_Schar2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
    p_rnode%dmemBytes = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_CHAR2BYTES)
  CASE DEFAULT
    PRINT *,'Error: unknown mem type'
    CALL sys_halt()
  END SELECT

  p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
  IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
    p_rheap%dtotalMemMax = p_rheap%dtotalMem

  ! Initialise the storage block
  CALL storage_initialiseBlock (ihandle, cinitNewBlock, rheap, Ilbound(2))

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_free (ihandle, rheap)

!<description>
  ! This routine releases a handle from a heap and deallocates the
  ! associated memory. ihandle is set to ST_NOHANDLE upon return.
!</description>

!<inputoutput>

  ! Handle of the memory block to be releases
  INTEGER :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .LE. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_free: Releasing ST_NOHANDLE is not allowed!'
    CALL sys_halt()
  END IF

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Is the node associated at all?
  IF (p_rnode%idataType .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_free: Trying to release nonexistent handle!'
    PRINT *,'Handle number: ',ihandle
    CALL sys_halt()
  END IF

  ! Release the memory assigned to that handle.
  IF (ASSOCIATED(p_rnode%p_Fsingle1D))  DEALLOCATE(p_rnode%p_Fsingle1D)
  IF (ASSOCIATED(p_rnode%p_Ddouble1D))  DEALLOCATE(p_rnode%p_Ddouble1D)
  IF (ASSOCIATED(p_rnode%p_Iinteger1D)) DEALLOCATE(p_rnode%p_Iinteger1D)
  IF (ASSOCIATED(p_rnode%p_Blogical1D)) DEALLOCATE(p_rnode%p_Blogical1D)
  IF (ASSOCIATED(p_rnode%p_Schar1D))    DEALLOCATE(p_rnode%p_Schar1D)
  IF (ASSOCIATED(p_rnode%p_Fsingle2D))  DEALLOCATE(p_rnode%p_Fsingle2D)
  IF (ASSOCIATED(p_rnode%p_Ddouble2D))  DEALLOCATE(p_rnode%p_Ddouble2D)
  IF (ASSOCIATED(p_rnode%p_Iinteger2D)) DEALLOCATE(p_rnode%p_Iinteger2D)
  IF (ASSOCIATED(p_rnode%p_Blogical2D)) DEALLOCATE(p_rnode%p_Blogical2D)
  IF (ASSOCIATED(p_rnode%p_Schar2D))    DEALLOCATE(p_rnode%p_Schar2D)

  ! Release the handle itself.
  CALL storage_releasehandle (ihandle,p_rheap)

  ! And finally reset the handle to ST_NOHANDLE.
  ihandle = ST_NOHANDLE

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_clear (ihandle, rheap)

!<description>
  ! This routine clears an array identified by ihandle; all entries are
  ! overwritten by 0.
!</description>

!<inputoutput>

  ! Handle of the memory block to be releases
  INTEGER :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .LE. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_clear: Handle invalid!'
    CALL sys_halt()
  END IF

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Is the node associated at all?
  IF (p_rnode%idataType .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_clear: Trying to release nonexistent handle!'
    PRINT *,'Handle number: ',ihandle
    CALL sys_halt()
  END IF

  ! What are we?
  SELECT CASE (p_rnode%idimension)
  CASE (1)
    SELECT CASE (p_rnode%idataType)
    CASE (ST_SINGLE)
      p_rnode%p_Fsingle1D = 0.0_SP
    CASE (ST_DOUBLE)
      p_rnode%p_Ddouble1D = 0.0_DP
    CASE (ST_INT)
      p_rnode%p_Iinteger1D = 0_I32
    CASE (ST_LOGICAL)
      p_rnode%p_Blogical1D = .FALSE.
    CASE (ST_CHAR)
      p_rnode%p_Schar1D = ACHAR(0)
    END SELECT
  CASE (2)
    SELECT CASE (p_rnode%idataType)
    CASE (ST_SINGLE)
      p_rnode%p_Fsingle2D = 0.0_SP
    CASE (ST_DOUBLE)
      p_rnode%p_Ddouble2D = 0.0_DP
    CASE (ST_INT)
      p_rnode%p_Iinteger2D = 0_I32
    CASE (ST_LOGICAL)
      p_rnode%p_Blogical2D = .FALSE.
    CASE (ST_CHAR)
      p_rnode%p_Schar2D = ACHAR(0)
    END SELECT
  CASE DEFAULT
    PRINT *,'storage_clear: invalid dimension.'
    CALL sys_halt()
  END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getsize1D (ihandle, isize, rheap)

!<description>
  ! Returns the length of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  INTEGER, INTENT(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!<output>
  ! Length of the array identified by ihandle.
  INTEGER(I32), INTENT(OUT) :: isize
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .LE. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getsize1D: Handle invalid!'
    CALL sys_halt()
  END IF

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Is the node associated at all?
  IF (p_rnode%idataType .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getsize1D: Handle invalid!'
    PRINT *,'Handle number: ',ihandle
    CALL sys_halt()
  END IF

  ! What are we?
  IF (p_rnode%idimension .NE. 1) THEN
    PRINT *,'Error in storage_getsize1D: Handle ',ihandle,' is not 1-dimensional!'
    CALL sys_halt()
  END IF

  SELECT CASE (p_rnode%idataType)
  CASE (ST_SINGLE)
    isize = SIZE(p_rnode%p_Fsingle1D)
  CASE (ST_DOUBLE)
    isize = SIZE(p_rnode%p_Ddouble1D)
  CASE (ST_INT)
    isize = SIZE(p_rnode%p_Iinteger1D)
  CASE (ST_LOGICAL)
    isize = SIZE(p_rnode%p_Blogical1D)
  CASE (ST_CHAR)
    isize = SIZE(p_rnode%p_Schar1D)
  CASE DEFAULT
    PRINT *,'Error in storage_getsize1D: Invalid data type!'
    CALL sys_halt()
  END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getsize2D (ihandle, isize, rheap)

!<description>
  ! Returns the length of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  INTEGER, INTENT(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!<output>
  ! Length of each dimension of the array identified by ihandle.
  INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Isize
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .LE. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getsize2D: Handle invalid!'
    CALL sys_halt()
  END IF

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Is the node associated at all?
  IF (p_rnode%idataType .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getsize2D: Handle invalid!'
    PRINT *,'Handle number: ',ihandle
    CALL sys_halt()
  END IF

  ! What are we?
  IF (p_rnode%idimension .NE. 2) THEN
    PRINT *,'Error in storage_getsize1D: Handle ',ihandle,' is not 2-dimensional!'
    CALL sys_halt()
  END IF

  SELECT CASE (p_rnode%idataType)
  CASE (ST_SINGLE)
    Isize = SHAPE(p_rnode%p_Fsingle2D)
  CASE (ST_DOUBLE)
    Isize = SHAPE(p_rnode%p_Ddouble2D)
  CASE (ST_INT)
    Isize = SHAPE(p_rnode%p_Iinteger2D)
  CASE (ST_LOGICAL)
    Isize = SHAPE(p_rnode%p_Blogical2D)
  CASE (ST_CHAR)
    Isize = SHAPE(p_rnode%p_Schar2D)
  CASE DEFAULT
    PRINT *,'Error in storage_getsize2D: Invalid data type!'
    CALL sys_halt()
  END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_int (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_int: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_INT) THEN
    PRINT *,'storage_getbase_int: Wrong data format!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_intUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_int: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_INT) THEN
    PRINT *,'storage_getbase_int: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D))) THEN
    PRINT *,'storage_getbase_intUBnd: Upper bound invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer
  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D(:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_intLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle
  
  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_int: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_INT) THEN
    PRINT *,'storage_getbase_int: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((lbnd .LT. 1) .OR. &
      (lbnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D)) .OR. &
      (ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D)) .OR. &
      (lbnd .GT. ubnd)) THEN
    PRINT *,'storage_getbase_intLUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D(lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_single (ihandle, p_Sarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(SP), DIMENSION(:), POINTER :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_single: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_SINGLE) THEN
    PRINT *,'storage_getbase_single: Wrong data format!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_singleUBnd (ihandle, p_Sarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(SP), DIMENSION(:), POINTER :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_single: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_SINGLE) THEN
    PRINT *,'storage_getbase_single: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D))) THEN
    PRINT *,'storage_getbase_singleUBnd: Upper bound invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D(:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_singleLUBnd (ihandle, p_Sarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(SP), DIMENSION(:), POINTER :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_single: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_SINGLE) THEN
    PRINT *,'storage_getbase_single: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((lbnd .LT. 1) .OR. &
      (lbnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D)) .OR. &
      (ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D)) .OR. &
      (lbnd .GT. ubnd)) THEN
    PRINT *,'storage_getbase_singleLUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D(lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_double (ihandle, p_Darray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(DP), DIMENSION(:), POINTER :: p_Darray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_double: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_DOUBLE) THEN

    PRINT *,'storage_getbase_double: Wrong data format!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_doubleUBnd (ihandle, p_Darray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(DP), DIMENSION(:), POINTER :: p_Darray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_double: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_DOUBLE) THEN
    PRINT *,'storage_getbase_double: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D))) THEN
    PRINT *,'storage_getbase_doubleUBnd: Upper bound invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer
  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D(:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_doubleLUBnd (ihandle, p_Darray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(DP), DIMENSION(:), POINTER :: p_Darray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_double: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_DOUBLE) THEN
    PRINT *,'storage_getbase_double: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((lbnd .LT. 1) .OR. &
      (lbnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D)) .OR. &
      (ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D)) .OR. &
      (lbnd .GT. ubnd)) THEN
    PRINT *,'storage_getbase_doubleLUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer
  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D(lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_logical (ihandle, p_Larray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  LOGICAL, DIMENSION(:), POINTER :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_logical: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_LOGICAL) THEN
    PRINT *,'storage_getbase_logical: Wrong data format!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_logicalUBnd (ihandle, p_Larray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  LOGICAL, DIMENSION(:), POINTER :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_logical: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_LOGICAL) THEN
    PRINT *,'storage_getbase_logical: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D))) THEN
    PRINT *,'storage_getbase_logicalUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D(:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_logicalLUBnd (ihandle, p_Larray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopts the given lower and upper bounds.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  LOGICAL, DIMENSION(:), POINTER :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_logical: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_LOGICAL) THEN
    PRINT *,'storage_getbase_logical: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((lbnd .LT. 1) .OR. &
      (lbnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .OR. &
      (ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .OR. &
      (lbnd .GT. ubnd)) THEN
    PRINT *,'storage_getbase_logicalUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D(lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_char (ihandle, p_Carray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  CHARACTER, DIMENSION(:), POINTER :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_char: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_CHAR) THEN
    PRINT *,'storage_getbase_char: Wrong data format!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_charUBnd (ihandle, p_Carray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle
  
  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  CHARACTER, DIMENSION(:), POINTER :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_char: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_CHAR) THEN
    PRINT *,'storage_getbase_char: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D))) THEN
    PRINT *,'storage_getbase_charLUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D(:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_charLUBnd (ihandle, p_Carray, lbnd, ubnd,rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopts the given lower and upper bounds.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  CHARACTER, DIMENSION(:), POINTER :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_char: Wrong handle'
    CALL sys_halt()
  END IF

  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_CHAR) THEN
    PRINT *,'storage_getbase_char: Wrong data format!'
    CALL sys_halt()
  END IF

  IF ((lbnd .LT. 1) .OR. &
      (lbnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .OR. &
      (ubnd .LT. 1) .OR. &
      (ubnd .GT. SIZE(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .OR. &
      (lbnd .GT. ubnd)) THEN
    PRINT *,'storage_getbase_charLUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D(lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_int2D (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_int2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_int2DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array and adopts the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_int2D: Wrong handle'
    CALL sys_halt()
  END IF

  IF ((ubnd .LT. 1) .OR. &
      (ubnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2))) THEN
    PRINT *,'storage_getbase_int2DUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D(:,:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_int2DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array and adopts the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_int2D: Wrong handle'
    CALL sys_halt()
  END IF

  IF ((lbnd .LT. 1) .OR. &
      (lbnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2)) .OR. &
      (ubnd .LT. 1) .OR. &
      (ubnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2)) .OR. &
      (lbnd .GT. ubnd)) THEN
    PRINT *,'storage_getbase_int2DLUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D(:,lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_single2D (ihandle, p_Sarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! single precision array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(SP), DIMENSION(:,:), POINTER :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_single2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_single2DUBnd (ihandle, p_Sarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single precision array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(SP), DIMENSION(:,:), POINTER :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_single2D: Wrong handle'
    CALL sys_halt()
  END IF

  IF ((ubnd .LT. 1) .OR. &
      (ubnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2))) THEN
    PRINT *,'storage_getbase_single2DUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D(:,:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_single2DLUBnd (ihandle, p_Sarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! single precision array and adopt the given lower and upper bounds
  ! for the second dimension.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd
  
  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  REAL(SP), DIMENSION(:,:), POINTER :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_single2D: Wrong handle'
    CALL sys_halt()
  END IF

  IF ((lbnd .LT. 1) .OR. &
      (lbnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2)) .OR. &
      (ubnd .LT. 1) .OR. &
      (ubnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2)) .OR. &
      (lbnd .GT. ubnd)) THEN
    PRINT *,'storage_getbase_single2DLUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D(:,lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_double2D (ihandle, p_Darray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array.

!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  REAL(DP), DIMENSION(:,:), POINTER :: p_Darray
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_double2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_double2DUBnd (ihandle, p_Darray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array and adopt the given upper bound for the second dimension.

!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  REAL(DP), DIMENSION(:,:), POINTER :: p_Darray
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_double2D: Wrong handle'
    CALL sys_halt()
  END IF

  IF ((ubnd .LT. 1) .OR. &
      (ubnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2))) THEN
    PRINT *,'storage_getbase_double2DUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D(:,:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_double2DLUBnd (ihandle, p_Darray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array and adopt the given lower and upper bounds
  ! for the second dimension.

!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  REAL(DP), DIMENSION(:,:), POINTER :: p_Darray
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_double2D: Wrong handle'
    CALL sys_halt()
  END IF

  IF ((lbnd .LT. 1) .OR. &
      (lbnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2)) .OR. &
      (ubnd .LT. 1) .OR. &
      (ubnd .GT. UBOUND(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2)) .OR. &
      (lbnd .GT. ubnd)) THEN
    PRINT *,'storage_getbase_double2DLUBnd: Bounds invalid invalid!'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D(:,lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_logical2D (ihandle, p_Larray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  LOGICAL, DIMENSION(:,:), POINTER :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_logical2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_logical2DUBnd (ihandle, p_Larray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  LOGICAL, DIMENSION(:,:), POINTER :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_logical2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D(:,:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_logical2DLUBnd (ihandle, p_Larray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopt the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  LOGICAL, DIMENSION(:,:), POINTER :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_logical2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D(:,lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_char2D (ihandle, p_Carray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  CHARACTER, DIMENSION(:,:), POINTER :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_char2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_char2DUBnd (ihandle, p_Carray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  CHARACTER, DIMENSION(:,:), POINTER :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_char2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D(:,:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_char2DLUBnd (ihandle, p_Carray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopt the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! The lower bound
  INTEGER, INTENT(IN) :: lbnd

  ! The upper bound
  INTEGER, INTENT(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  CHARACTER, DIMENSION(:,:), POINTER :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .EQ. ST_NOHANDLE) THEN
    PRINT *,'storage_getbase_char2D: Wrong handle'
    CALL sys_halt()
  END IF

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D(:,lbnd:ubnd)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_copy(h_source, h_dest, rheap)

!<description>
  ! This routine copies the information of one array to another.
  ! The structure of the arrays behind h_source and h_dest must be the
  ! same!
!</description>

!<input>
  ! Handle of the source array to copy
  INTEGER, INTENT(IN) :: h_source

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(IN), TARGET, OPTIONAL :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  INTEGER, INTENT(INOUT) :: h_dest
!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rsource, p_rdest
  INTEGER(I32) :: i,j
  INTEGER(I32), DIMENSION(2) :: Isize

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbase
    END IF

    IF (h_source .EQ. ST_NOHANDLE) THEN
      PRINT *,'storage_copy: Wrong handle'
      CALL sys_halt()
    END IF
    IF (.NOT. ASSOCIATED(p_rheap%p_Rdescriptors)) THEN
      PRINT *,'storage_copy: Heap not initialised!'
      CALL sys_halt()
    END IF

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    IF (h_dest .EQ. ST_NOHANDLE) THEN
      ! Create a new array in the same size and structure
      ! as h_source.
      SELECT CASE (p_rsource%idimension)
      CASE (1)
        SELECT CASE (p_rsource%idataType)
        CASE (ST_DOUBLE)
          CALL storage_new ('storage_copy',p_rsource%sname,&
                            INT(SIZE(p_rsource%p_Ddouble1D),I32),&
                            ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_SINGLE)
          CALL storage_new ('storage_copy',p_rsource%sname,&
                            INT(SIZE(p_rsource%p_Fsingle1D),I32),&
                            ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_INT)
          CALL storage_new ('storage_copy',p_rsource%sname,&
                            INT(SIZE(p_rsource%p_Iinteger1D),I32),&
                            ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_LOGICAL)
          CALL storage_new ('storage_copy',p_rsource%sname,&
                            INT(SIZE(p_rsource%p_Blogical1D),I32),&
                            ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_CHAR)
          CALL storage_new ('storage_copy',p_rsource%sname,&
                            INT(SIZE(p_rsource%p_Schar1D),I32),&
                            ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        END SELECT
      CASE (2)
        SELECT CASE (p_rsource%IdataType)
        CASE (ST_DOUBLE)
          Isize = UBOUND(p_rsource%p_Ddouble2D)   ! =SIZE(...) here
          CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_SINGLE)
          Isize = UBOUND(p_rsource%p_Fsingle2D)   ! =SIZE(...) here
          CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_INT)
          Isize = UBOUND(p_rsource%p_Iinteger2D)      ! =SIZE(...) here
          CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_LOGICAL)
          Isize = UBOUND(p_rsource%p_Blogical2D)      ! =SIZE(...) here
          CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_CHAR)
          Isize = UBOUND(p_rsource%p_Schar2D)      ! =SIZE(...) here
          CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        END SELECT
      END SELECT

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it's correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      
    END IF

    p_rdest => p_rheap%p_Rdescriptors(h_dest)

    ! 1D/2D the same?
    IF (p_rsource%idimension .NE. p_rdest%idimension) THEN
      PRINT *,'storage_copy: Structure different!'
      CALL sys_halt()
    END IF

    ! What is to copy
    SELECT CASE (p_rsource%idimension)
    CASE (1)
      SELECT CASE (p_rsource%idataType)
      CASE (ST_DOUBLE)
        SELECT CASE (p_rdest%idataType)
        CASE (ST_DOUBLE)
          CALL lalg_copyVectorDble (p_rsource%p_Ddouble1D,p_rdest%p_Ddouble1D)
        CASE (ST_SINGLE)
          CALL lalg_copyVectorDblSngl (p_rsource%p_Ddouble1D,p_rdest%p_Fsingle1D)
        CASE DEFAULT
          PRINT *,'storage_copy: Unsupported data type combination'
          CALL sys_halt()
        END SELECT

      CASE (ST_SINGLE)
        SELECT CASE (p_rdest%idataType)
        CASE (ST_DOUBLE)
          CALL lalg_copyVectorSnglDbl (p_rsource%p_Fsingle1D,p_rdest%p_Ddouble1D)
        CASE (ST_SINGLE)
          CALL lalg_copyVectorSngl (p_rsource%p_Fsingle1D,p_rdest%p_Fsingle1D)
        CASE DEFAULT
          PRINT *,'storage_copy: Unsupported data type combination'
          CALL sys_halt()
        END SELECT

      CASE (ST_INT)
        IF (p_rdest%idataType .EQ. ST_INT) THEN
          CALL lalg_copyVectorInt (p_rsource%p_Iinteger1D,p_rdest%p_Iinteger1D)
        ELSE
          PRINT *,'storage_copy: Unsupported data type combination'
          CALL sys_halt()
        END IF

      CASE (ST_LOGICAL)
        IF (p_rdest%idataType .EQ. ST_LOGICAL) THEN
          ! Copy by hand
          DO i=1,SIZE(p_rsource%p_Blogical1D)
            p_rdest%p_Blogical1D(i) = p_rsource%p_Blogical1D(i)
          END DO
        ELSE
          PRINT *,'storage_copy: Unsupported data type combination'
          CALL sys_halt()
        END IF

      CASE (ST_CHAR)
        IF (p_rdest%idataType .EQ. ST_CHAR) THEN
          ! Copy by hand
          DO i=1,SIZE(p_rsource%p_Schar1D)
            p_rdest%p_Schar1D(i) = p_rsource%p_Schar1D(i)
          END DO
        ELSE
          PRINT *,'storage_copy: Unsupported data type combination'
          CALL sys_halt()
        END IF

      CASE DEFAULT
        PRINT *,'storage_copy: Unknown data type'
        CALL sys_halt()
      END SELECT

    CASE (2)
      SELECT CASE (p_rsource%idataType)
      CASE (ST_DOUBLE)
        IF ((UBOUND(p_rsource%p_Ddouble2D,1) .NE. UBOUND(p_rdest%p_Ddouble2D,1)) .OR.&
            (UBOUND(p_rsource%p_Ddouble2D,2) .NE. UBOUND(p_rdest%p_Ddouble2D,2))) THEN
          PRINT *,'storage_copy: Structure different!'
          CALL sys_halt()
        END IF

        SELECT CASE (p_rdest%idataType)
        CASE (ST_DOUBLE)
          ! Copy by hand
          DO j=1,UBOUND(p_rsource%p_Ddouble2D,2)
            DO i=1,UBOUND(p_rsource%p_Ddouble2D,1)
              p_rdest%p_Ddouble2D(i,j) = p_rsource%p_Ddouble2D(i,j)
            END DO
          END DO

        CASE (ST_SINGLE)
          ! Copy by hand
          DO j=1,UBOUND(p_rsource%p_Fsingle2D,2)
            DO i=1,UBOUND(p_rsource%p_Fsingle2D,1)
              p_rdest%p_Fsingle2D(i,j) = p_rsource%p_Ddouble2D(i,j)
            END DO
          END DO

        CASE (ST_INT)
          ! Copy by hand
          DO j=1,UBOUND(p_rsource%p_Iinteger2D,2)
            DO i=1,UBOUND(p_rsource%p_Iinteger2D,1)
              p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Ddouble2D(i,j)
            END DO
          END DO

        CASE DEFAULT ! Might be ST_LOGICAL or ST_CHAR
          PRINT *,'storage_copy: Unsupported data type combination'
          CALL sys_halt()

        END SELECT

      CASE (ST_SINGLE)
        IF ((UBOUND(p_rsource%p_Fsingle2D,1) .NE. UBOUND(p_rdest%p_Fsingle2D,1)) .OR.&
            (UBOUND(p_rsource%p_Fsingle2D,2) .NE. UBOUND(p_rdest%p_Fsingle2D,2))) THEN
          PRINT *,'storage_copy: Structure different!'
          CALL sys_halt()
        END IF

        SELECT CASE (p_rdest%idataType)
        CASE (ST_DOUBLE)
          ! Copy by hand
          DO j=1,UBOUND(p_rsource%p_Ddouble2D,2)
            DO i=1,UBOUND(p_rsource%p_Ddouble2D,1)
              p_rdest%p_Ddouble2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            END DO
          END DO

        CASE (ST_SINGLE)
          ! Copy by hand
          DO j=1,UBOUND(p_rsource%p_Fsingle2D,2)
            DO i=1,UBOUND(p_rsource%p_Fsingle2D,1)
              p_rdest%p_Fsingle2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            END DO
          END DO

        CASE (ST_INT)
          ! Copy by hand
          DO j=1,UBOUND(p_rsource%p_Iinteger2D,2)
            DO i=1,UBOUND(p_rsource%p_Iinteger2D,1)
              p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            END DO
          END DO

        ! Might be ST_LOGICAL or ST_CHAR
        CASE DEFAULT
          PRINT *,'storage_copy: Unsupported data type combination'
          CALL sys_halt()

        END SELECT

      CASE (ST_INT)
        IF ((UBOUND(p_rsource%p_Iinteger2D,1) .NE. UBOUND(p_rdest%p_Iinteger2D,1)) .OR.&
            (UBOUND(p_rsource%p_Iinteger2D,2) .NE. UBOUND(p_rdest%p_Iinteger2D,2))) THEN
          PRINT *,'storage_copy: Structure different!'
          CALL sys_halt()
        END IF
        IF (p_rdest%idataType .NE. ST_INT) THEN
          PRINT *,'storage_copy: unsupported data type combination'
          CALL sys_halt()
        END IF

        ! Copy by hand
        DO j=1,UBOUND(p_rsource%p_Iinteger2D,2)
          DO i=1,UBOUND(p_rsource%p_Iinteger2D,1)
            p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Iinteger2D(i,j)
          END DO
        END DO

      CASE (ST_LOGICAL)
        IF ((UBOUND(p_rsource%p_Blogical2D,1) .NE. UBOUND(p_rdest%p_Blogical2D,1)) .OR.&
            (UBOUND(p_rsource%p_Blogical2D,2) .NE. UBOUND(p_rdest%p_Blogical2D,2))) THEN
          PRINT *,'storage_copy: Structure different!'
          CALL sys_halt()
        END IF
        IF (p_rdest%idataType .NE. ST_LOGICAL) THEN
          PRINT *,'storage_copy: unsupported data type combination'
          CALL sys_halt()
        END IF

        ! Copy by hand
        DO j=1,UBOUND(p_rsource%p_Blogical2D,2)
          DO i=1,UBOUND(p_rsource%p_Blogical2D,1)
            p_rdest%p_Blogical2D(i,j) = p_rsource%p_Blogical2D(i,j)
          END DO
        END DO

      CASE (ST_CHAR)
        IF ((UBOUND(p_rsource%p_Schar2D,1) .NE. UBOUND(p_rdest%p_Schar2D,1)) .OR.&
            (UBOUND(p_rsource%p_Schar2D,2) .NE. UBOUND(p_rdest%p_Schar2D,2))) THEN
          PRINT *,'storage_copy: Structure different!'
          CALL sys_halt()
        END IF
        IF (p_rdest%idataType .NE. ST_CHAR) THEN
          PRINT *,'storage_copy: unsupported data type combination'
          CALL sys_halt()
        END IF

        ! Copy by hand
        DO j=1,UBOUND(p_rsource%p_Schar2D,2)
          DO i=1,UBOUND(p_rsource%p_Schar2D,1)
            p_rdest%p_Schar2D(i,j) = p_rsource%p_Schar2D(i,j)
          END DO
        END DO

      CASE DEFAULT
        PRINT *,'storage_copy: Unknown data type'
        CALL sys_halt()
      END SELECT
    END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_copy_explicit(h_source, h_dest, istart_source, &
             istart_dest, ilength, rheap)

!<description>
  ! This routine copies the information of one array to another.
  ! The structure of the arrays behind h_source and h_dest must be
  ! "similar" in the following sense: Both datatypes and dimensions
  ! must be the same. The routine allows for copying only parts of
  ! of the arrays. Therefor the relevant parts of the arrays must
  ! be the same!
!</description>

!<input>
  ! Handle of the source array to copy
  INTEGER, INTENT(IN) :: h_source

  ! First entry of the source array to copy
  INTEGER, INTENT(IN) :: istart_source

  ! First entry of the destination array where to copy
  INTEGER, INTENT(IN) :: istart_dest

  ! Length of the array to copy
  INTEGER, INTENT(IN) :: ilength

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(IN), TARGET, OPTIONAL :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  INTEGER, INTENT(INOUT) :: h_dest
!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    TYPE(t_storageBlock), POINTER :: p_rheap
    TYPE(t_storageNode), POINTER :: p_rsource, p_rdest
    INTEGER(I32) :: i
    CHARACTER(LEN=SYS_NAMELEN) :: sname = ''

    ! Check if the start address is positive
    IF (istart_source <= 0 .OR. istart_dest <= 0) THEN
      PRINT *, 'storage_copy_explicit: start address must be positive'
      CALL sys_halt()
    END IF

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbase
    END IF

    IF (h_source .EQ. ST_NOHANDLE) THEN
      PRINT *,'storage_copy_explicit: Wrong handle'
      CALL sys_halt()
    END IF
    IF (.NOT. ASSOCIATED(p_rheap%p_Rdescriptors)) THEN
      PRINT *,'storage_copy_explicit: Heap not initialised!'
      CALL sys_halt()
    END IF

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    IF (h_dest .EQ. ST_NOHANDLE) THEN
      ! Create a new array in the same size and structure
      ! as h_source.
      IF (p_rsource%idimension /= 1) THEN
        PRINT *, 'storage_copy_explicit: only 1D arrays are allowed'
        CALL sys_halt()
      END IF

      SELECT CASE (p_rsource%idataType)
      CASE (ST_DOUBLE)
         CALL storage_new ('storage_copy_explicit',p_rsource%sname,&
              INT(SIZE(p_rsource%p_Ddouble1D),I32),&
              ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      CASE (ST_SINGLE)
         CALL storage_new ('storage_copy_explicit',p_rsource%sname,&
              INT(SIZE(p_rsource%p_Fsingle1D),I32),&
              ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      CASE (ST_INT)
         CALL storage_new ('storage_copy_explicit',p_rsource%sname,&
              INT(SIZE(p_rsource%p_Iinteger1D),I32),&
              ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      CASE (ST_LOGICAL)
         CALL storage_new ('storage_copy_explicit',p_rsource%sname,&
              INT(SIZE(p_rsource%p_Blogical1D),I32),&
              ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      CASE (ST_CHAR)
         CALL storage_new ('storage_copy_explicit',p_rsource%sname,&
              INT(SIZE(p_rsource%p_Schar1D),I32),&
              ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      END SELECT

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it's correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      
    END IF

    p_rdest => p_rheap%p_Rdescriptors(h_dest)

    ! 1D/2D the same?
    IF (p_rsource%idimension .NE. p_rdest%idimension) THEN
      PRINT *,'storage_copy_explicit: Structure different!'
      CALL sys_halt()
    END IF

    ! What is to copy
    SELECT CASE (p_rsource%idataType)
    CASE (ST_DOUBLE)
       SELECT CASE (p_rdest%idataType)
       CASE (ST_DOUBLE)
          IF (istart_source+ilength-1 > SIZE(p_rsource%p_Ddouble1D) .OR. &
               istart_dest+ilength-1 > SIZE(p_rdest%p_Ddouble1D)) THEN
            PRINT *, 'storage_copy_explicit: Subarrays incompatible!'
            CALL sys_halt()
          END IF
          ! Copy by hand
          DO i=1,ilength
            p_rdest%p_Ddouble1D(istart_dest+i-1) = &
              p_rsource%p_Ddouble1D(istart_source+i-1)
          END DO
       CASE (ST_SINGLE)
          IF (istart_source+ilength-1 > SIZE(p_rsource%p_Ddouble1D) .OR. &
               istart_dest+ilength-1 > SIZE(p_rdest%p_Fsingle1D)) THEN
            PRINT *, 'storage_copy_explicit: Subarrays incompatible!'
            CALL sys_halt()
          END IF
          ! Copy by hand
          DO i=1,ilength
            p_rdest%p_Fsingle1D(istart_dest+i-1) = &
              p_rsource%p_Ddouble1D(istart_source+i-1)
          END DO
       CASE DEFAULT
          PRINT *,'storage_copy_explicit: Unsupported data type combination'
          CALL sys_halt()
       END SELECT

    CASE (ST_SINGLE)
       SELECT CASE (p_rdest%idataType)
       CASE (ST_DOUBLE)
          IF (istart_source+ilength-1 > SIZE(p_rsource%p_Fsingle1D) .OR. &
               istart_dest+ilength-1 > SIZE(p_rdest%p_Ddouble1D)) THEN
             PRINT *, 'storage_copy_explicit: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO i=1,ilength
             p_rdest%p_Ddouble1D(istart_dest+i-1) = &
               p_rsource%p_Fsingle1D(istart_source+i-1)
          END DO
       CASE (ST_SINGLE)
          IF (istart_source+ilength-1 > SIZE(p_rsource%p_Fsingle1D) .OR. &
               istart_dest+ilength-1 > SIZE(p_rdest%p_Fsingle1D)) THEN
             PRINT *, 'storage_copy_explicit: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO i=1,ilength
             p_rdest%p_Fsingle1D(istart_dest+i-1) = &
               p_rsource%p_Fsingle1D(istart_source+i-1)
          END DO
       CASE DEFAULT
          PRINT *,'storage_copy_explicit: Unsupported data type combination'
          CALL sys_halt()
       END SELECT

    CASE (ST_INT)
       IF (p_rdest%idataType .EQ. ST_INT) THEN
          IF (istart_source+ilength-1 > SIZE(p_rsource%p_Iinteger1D) .OR. &
               istart_dest+ilength-1 > SIZE(p_rdest%p_Iinteger1D)) THEN
             PRINT *, 'storage_copy_explicit: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO i=1,ilength
             p_rdest%p_Iinteger1D(istart_dest+i-1) = p_rsource%p_Iinteger1D(istart_source+i-1)
          END DO
       ELSE
          PRINT *,'storage_copy_explicit: Unsupported data type combination'
          CALL sys_halt()
       END IF

    CASE (ST_LOGICAL)
       IF (p_rdest%idataType .EQ. ST_LOGICAL) THEN
          IF (istart_source+ilength-1 > SIZE(p_rsource%p_Blogical1D) .OR. &
               istart_dest+ilength-1 > SIZE(p_rdest%p_Blogical1D)) THEN
             PRINT *, 'storage_copy_explicit: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO i=1,ilength
             p_rdest%p_Blogical1D(istart_dest+i-1) = p_rsource%p_Blogical1D(istart_source+i-1)
          END DO
       ELSE
          PRINT *,'storage_copy_explicit: Unsupported data type combination'
          CALL sys_halt()
       END IF

    CASE (ST_CHAR)
       IF (p_rdest%idataType .EQ. ST_CHAR) THEN
          IF (istart_source+ilength-1 > SIZE(p_rsource%p_Schar1D) .OR. &
               istart_dest+ilength-1 > SIZE(p_rdest%p_Schar1D)) THEN
             PRINT *, 'storage_copy_explicit: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO i=1,ilength
             p_rdest%p_Schar1D(istart_dest+i-1) = p_rsource%p_Schar1D(istart_source+i-1)
          END DO
       ELSE
          PRINT *,'storage_copy_explicit: Unsupported data type combination'
          CALL sys_halt()
       END IF

    CASE DEFAULT
       PRINT *,'storage_copy_explicit: Unknown data type'
       CALL sys_halt()
    END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_copy_explicit2D(h_source, h_dest, Istart_source, &
      Istart_dest, Ilength, rheap)

!<description>
  ! This routine copies the information of one array to another.
  ! The structure of the arrays behind h_source and h_dest must be
  ! "similar" in the following sense: Both datatypes and dimensions
  ! must be the same. The routine allows for copying only parts of
  ! of the arrays. Therefor the relevant parts of the arrays must
  ! be the same!
!</description>

!<input>
  ! Handle of the source array to copy
  INTEGER, INTENT(IN) :: h_source

  ! First entry of the source array to copy
  INTEGER, DIMENSION(2), INTENT(IN) :: Istart_source

  ! First entry of the destination array where to copy
  INTEGER, DIMENSION(2), INTENT(IN) :: Istart_dest

  ! Length of the array to copy
  INTEGER, DIMENSION(2), INTENT(IN) :: Ilength

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(IN), TARGET, OPTIONAL :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  INTEGER, INTENT(INOUT) :: h_dest
!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rsource, p_rdest
  INTEGER(I32) :: i,j
  INTEGER(I32), DIMENSION(2) :: Isize

    ! Check if the start address is positive
    IF (ANY(istart_source <= 0) .OR. ANY(istart_dest <= 0)) THEN
      PRINT *, 'storage_copy_explicit: start address must be positive'
      CALL sys_halt()
    END IF

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbase
    END IF

    IF (h_source .EQ. ST_NOHANDLE) THEN
      PRINT *,'storage_copy_explicit: Wrong handle'
      CALL sys_halt()
    END IF
    IF (.NOT. ASSOCIATED(p_rheap%p_Rdescriptors)) THEN
      PRINT *,'storage_copy_explicit: Heap not initialised!'
      CALL sys_halt()
    END IF

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    IF (h_dest .EQ. ST_NOHANDLE) THEN
      ! Create a new array in the same size and structure
      ! as h_source.
      IF (p_rsource%idimension /= 2) THEN
        PRINT *, 'storage_copy_explicit: only 1D arrays are allowed'
        CALL sys_halt()
      END IF

      SELECT CASE (p_rsource%IdataType)
      CASE (ST_DOUBLE)
         Isize = UBOUND(p_rsource%p_Ddouble2D)   ! =SIZE(...) here
         CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      CASE (ST_SINGLE)
         Isize = UBOUND(p_rsource%p_Fsingle2D)   ! =SIZE(...) here
         CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      CASE (ST_INT)
         Isize = UBOUND(p_rsource%p_Iinteger2D)      ! =SIZE(...) here
         CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      CASE (ST_LOGICAL)
         Isize = UBOUND(p_rsource%p_Blogical2D)      ! =SIZE(...) here
         CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      CASE (ST_CHAR)
         Isize = UBOUND(p_rsource%p_Schar2D)      ! =SIZE(...) here
         CALL storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      END SELECT

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it's correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      
    END IF

    p_rdest => p_rheap%p_Rdescriptors(h_dest)

    ! 1D/2D the same?
    IF (p_rsource%idimension .NE. p_rdest%idimension) THEN
      PRINT *,'storage_copy_explicit: Structure different!'
      CALL sys_halt()
    END IF

    ! What is to copy
    SELECT CASE (p_rsource%idataType)
    CASE (ST_DOUBLE)
       SELECT CASE (p_rdest%idataType)
       CASE (ST_DOUBLE)
          IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Ddouble2D,1) .OR. &
               Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Ddouble2D,2) .OR. &
               Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Ddouble2D,1) .OR. &
               Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Ddouble2D,2)) THEN
             PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO j=1,Ilength(2)
             DO i=1,Ilength(1)
                p_rdest%p_Ddouble2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             END DO
          END DO

       CASE (ST_SINGLE)
          IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Ddouble2D,1) .OR. &
               Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Ddouble2D,2) .OR. &
               Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Fsingle2D,1) .OR. &
               Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Fsingle2D,2)) THEN
             PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO j=1,Ilength(2)
             DO i=1,Ilength(1)
                p_rdest%p_Fsingle2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             END DO
          END DO

       CASE (ST_INT)
          IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Ddouble2D,1) .OR. &
               Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Ddouble2D,2) .OR. &
               Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Iinteger2D,1) .OR. &
               Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Iinteger2D,2)) THEN
             PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO j=1,Ilength(2)
             DO i=1,Ilength(1)
                p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             END DO
          END DO

       CASE DEFAULT ! Might be ST_LOGICAL or ST_CHAR
          PRINT *,'storage_copy_explicit2D: Unsupported data type combination'
          CALL sys_halt()

       END SELECT

    CASE (ST_SINGLE)
       SELECT CASE (p_rdest%idataType)
       CASE (ST_DOUBLE)
          IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Fsingle2D,1) .OR. &
               Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Fsingle2D,2) .OR. &
               Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Ddouble2D,1) .OR. &
               Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Ddouble2D,2)) THEN
             PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO j=1,Ilength(2)
             DO i=1,Ilength(1)
                p_rdest%p_Ddouble2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             END DO
          END DO

       CASE (ST_SINGLE)
          IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Fsingle2D,1) .OR. &
               Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Fsingle2D,2) .OR. &
               Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Fsingle2D,1) .OR. &
               Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Fsingle2D,2)) THEN
             PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO j=1,Ilength(2)
             DO i=1,Ilength(1)
                p_rdest%p_Fsingle2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             END DO
          END DO

       CASE (ST_INT)
          IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Fsingle2D,1) .OR. &
               Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Fsingle2D,2) .OR. &
               Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Iinteger2D,1) .OR. &
               Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Iinteger2D,2)) THEN
             PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
             CALL sys_halt()
          END IF
          ! Copy by hand
          DO j=1,Ilength(2)
             DO i=1,Ilength(1)
                p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             END DO
          END DO

       CASE DEFAULT ! Might be ST_LOGICAL or ST_CHAR
          PRINT *,'storage_copy_explicit2D: Unsupported data type combination'
          CALL sys_halt()

       END SELECT

    CASE (ST_INT)
       IF (p_rdest%idataType .NE. ST_INT) THEN
          PRINT *,'storage_copy_explicit2D: unsupported data type combination'
          CALL sys_halt()
       END IF
       IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Iinteger2D,1) .OR. &
            Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Iinteger2D,2) .OR. &
            Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Iinteger2D,1) .OR. &
            Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Iinteger2D,2)) THEN
          PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
          CALL sys_halt()
       END IF

       ! Copy by hand
       DO j=1,Ilength(2)
          DO i=1,Ilength(1)
             p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                  p_rsource%p_Iinteger2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          END DO
       END DO

    CASE (ST_LOGICAL)
       IF (p_rdest%idataType .NE. ST_LOGICAL) THEN
          PRINT *,'storage_copy_explicit2D: unsupported data type combination'
          CALL sys_halt()
       END IF
       IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Blogical2D,1) .OR. &
            Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Blogical2D,2) .OR. &
            Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Blogical2D,1) .OR. &
            Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Blogical2D,2)) THEN
          PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
          CALL sys_halt()
       END IF

       ! Copy by hand
       DO j=1,Ilength(2)
          DO i=1,Ilength(1)
             p_rdest%p_Blogical2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                  p_rsource%p_Blogical2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          END DO
       END DO

    CASE (ST_CHAR)
       IF (p_rdest%idataType .NE. ST_CHAR) THEN
          PRINT *,'storage_copy_explicit2D: unsupported data type combination'
          CALL sys_halt()
       END IF
       IF (Istart_source(1)+Ilength(1)-1 > UBOUND(p_rsource%p_Schar2D,1) .OR. &
            Istart_source(2)+Ilength(2)-1 > UBOUND(p_rsource%p_Schar2D,2) .OR. &
            Istart_dest(1)+Ilength(1)-1 > UBOUND(p_rdest%p_Schar2D,1) .OR. &
            Istart_dest(2)+Ilength(2)-1 > UBOUND(p_rdest%p_Schar2D,2)) THEN
          PRINT *, 'storage_copy_explicit2D: Subarrays incompatible!'
          CALL sys_halt()
       END IF

       ! Copy by hand
       DO j=1,Ilength(2)
          DO i=1,Ilength(1)
             p_rdest%p_Schar2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                  p_rsource%p_Schar2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          END DO
       END DO

    CASE DEFAULT
       PRINT *,'storage_copy: Unknown data type'
       CALL sys_halt()
    END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_info(bprintHandles,rheap)

!<description>
  ! This routine prints information about the current memory consumption
  ! in a memory block to screen.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the handles still remaining in the
  ! heap together with their names are printed to the terminal.
  LOGICAL, INTENT(IN), OPTIONAL :: bprintHandles

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(IN), TARGET, OPTIONAL :: rheap
!</input>

!</subroutine>

  ! local variables
  INTEGER :: i

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbase
    END IF

    CALL output_line ('Heap statistics:')
    CALL output_line ('----------------')
    IF (PRESENT(bprintHandles)) THEN
      IF (bprintHandles .AND. (p_rheap%ihandlesInUse .GT. 0)) THEN
        CALL output_line ('Handles on the heap: ')
        CALL output_lbrk ()
        ! Loop through the heap and search allocated handles
        DO i=1,SIZE(p_rheap%p_IfreeHandles)
          IF (p_rheap%p_Rdescriptors(i)%idataType .NE. ST_NOHANDLE) THEN
            IF (p_rheap%p_Rdescriptors(i)%idimension .EQ. 1) THEN
              CALL output_line ( &
                   'Handle ' // TRIM(sys_siL(i,10)) // ', 1D, Length=' // &
                   TRIM(sys_siL(INT(p_rheap%p_Rdescriptors(i)%dmemBytes),15)) //&
                   ', Type=' // TRIM(sys_siL(p_rheap%p_Rdescriptors(i)%idataType,15)) //&
                   ' Name=' // TRIM(ADJUSTL(p_rheap%p_Rdescriptors(i)%sname)) )
            ELSE
              CALL output_line ( &
                   'Handle ' // TRIM(sys_siL(i,10)) // ', 2D, Length=' // &
                   TRIM(sys_siL(INT(p_rheap%p_Rdescriptors(i)%dmemBytes),15)) // &
                   ', Type=' // TRIM(sys_siL(p_rheap%p_Rdescriptors(i)%idataType,15)) //&
                   ' Name=' // TRIM(ADJUSTL(p_rheap%p_Rdescriptors(i)%sname)) )
            END IF
          END IF
        END DO
        CALL output_lbrk ()
      END IF
    END IF

    CALL output_line ('Number of handles in use:        '//&
                      TRIM(sys_siL(p_rheap%ihandlesInUse,15)))
    IF (p_rheap%dtotalMem .GT. REAL(HUGE(0),DP)) THEN
      CALL output_line ('Memory in use (bytes):           '//&
                        TRIM(sys_sdL(p_rheap%dtotalMem,0)))
    ELSE
      CALL output_line ('Memory in use (bytes):           '//&
                        TRIM(sys_siL(INT(p_rheap%dtotalMem),15)))
    END IF
    CALL output_line ('Current total number of handles: '//&
                      TRIM(sys_siL(SIZE(p_rheap%p_IfreeHandles),15)))
    CALL output_line ('Maximum number of handles used:  '//&
                      TRIM(sys_siL(p_rheap%nhandlesInUseMax,15)))

    IF (p_rheap%dtotalMem .GT. REAL(HUGE(0),DP)) THEN
      CALL output_line ('Maximum used memory (bytes):     '//&
                        TRIM(sys_sdL(p_rheap%dtotalMemMax,0)))
    ELSE
      CALL output_line ('Maximum used memory (bytes):     '//&
                        TRIM(sys_siL(INT(p_rheap%dtotalMemMax),15)))
    END IF
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getdatatype (ihandle, idatatype, rheap)

!<description>
  ! Returns the datatype of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  INTEGER, INTENT(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!<output>
  ! Datatype of the array identified by ihandle.
  INTEGER(I32), INTENT(OUT) :: idatatype
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .LE. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getdatatype: Handle invalid!'
    CALL sys_halt()
  END IF

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  SELECT CASE (p_rnode%idataType)
  CASE (ST_SINGLE)
    idatatype = ST_SINGLE
  CASE (ST_DOUBLE)
    idatatype = ST_DOUBLE
  CASE (ST_INT)
    idatatype = ST_INT
  CASE (ST_LOGICAL)
    idatatype = ST_LOGICAL
  CASE (ST_CHAR)
    idatatype = ST_CHAR
  CASE (ST_NOHANDLE)
    PRINT *,'Error in storage_getdatatype: Handle invalid!'
    PRINT *,'Handle number: ',ihandle
    CALL sys_halt()
  CASE DEFAULT
    PRINT *,'Error in storage_getdatatype: Invalid data type!'
    CALL sys_halt()
  END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getdimension (ihandle, idimension, rheap)

!<description>
  ! Returns the dimension of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  INTEGER, INTENT(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!<output>
  ! Dimension of the array identified by ihandle.
  INTEGER(I32), INTENT(OUT) :: idimension
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode

  ! Get the heap to use - local or global one.

  IF(PRESENT(rheap)) THEN
    p_rheap => rheap
  ELSE
    p_rheap => rbase
  END IF

  IF (ihandle .LE. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getdatatype: Handle invalid!'
    CALL sys_halt()
  END IF

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  IF (ihandle .LE. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getsize1D: Handle invalid!'
    CALL sys_halt()
  END IF

  idimension = p_rnode%idimension

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_realloc (scall, isize, ihandle, cinitNewBlock, bcopy, rheap)

!<description>
  ! This routine reallocates an existing memory block wih a new desired
  ! size. In case of a multiple-dimension block, the last dimension
  ! is changed. isize is the size of the new memory block / the new size
  ! of the last dimension.
  !
  ! Warning: Reallocation of an array destroys all pointers associated with
  ! the corresponding handle!
!</description>

!<input>

  ! Name of the calling routine
  CHARACTER(LEN=*), INTENT(IN) :: scall

  ! Requested storage size for the memory block / the new size of the last
  ! dimension in the memory block identified by ihandle
  INTEGER(I32), INTENT(IN) :: isize

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO,
  ! ST_NEWBLOCK_NOINIT, ST_NEWBLOCK_ORDERED).
  ! Specifies how to initialise memory if isize > original array size.
  INTEGER, INTENT(IN) :: cinitNewBlock

  ! OPTIONAL: Copy old data.
  ! =TRUE: Copy data of old array to the new one.
  ! =FALSE: Reallocate memory, don't copy old data.
  ! If not specified, TRUE is assumed.
  LOGICAL, INTENT(IN), OPTIONAL :: bcopy

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

  ! Handle of the memory block.
  INTEGER, INTENT(INOUT) :: ihandle

!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    TYPE(t_storageBlock), POINTER :: p_rheap
    TYPE(t_storageNode), POINTER :: p_rnode

    ! New storage node
    TYPE(t_storageNode) :: rstorageNode

    ! size of the old 1-dimensional array
    INTEGER(I32) :: isizeOld

    ! size of the 1-dimensional array to be copied
    INTEGER :: isizeCopy

    ! size of the old 2-dimensional array
    INTEGER(I32), DIMENSION(2) :: Isize2Dold

    INTEGER(I32) :: i,j

    LOGICAL :: bcopyData

    IF (isize .EQ. 0) THEN
      ! Ok, not much to do...
      CALL storage_free(ihandle,rheap)
      RETURN
    END IF

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbase
    END IF

    ! Copy old data?

    IF (PRESENT(bcopy)) THEN
      bcopyData = bcopy
    ELSE
      bcopyData = .TRUE.
    END IF

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Copy the data of the old storage node to rstorageNode. That way
    ! we prepare a new storage node and will replace the old.
    rstorageNode = p_rnode

    ! Are we 1D or 2D?
    SELECT CASE(p_rnode%idimension)

    CASE (1)

      ! Get the size of the old storage node.
      SELECT CASE (p_rnode%idataType)
      CASE (ST_SINGLE)
        isizeOld = SIZE(p_rnode%p_Fsingle1D)
      CASE (ST_DOUBLE)
        isizeOld = SIZE(p_rnode%p_Ddouble1D)
      CASE (ST_INT)
        isizeOld = SIZE(p_rnode%p_Iinteger1D)
      CASE (ST_LOGICAL)
        isizeOld = SIZE(p_rnode%p_Blogical1D)
      CASE (ST_CHAR)
        isizeOld = SIZE(p_rnode%p_Schar1D)
      END SELECT

      ! Do we really have to change anything?
      IF (isize == isizeOld) RETURN

      ! Allocate new memory and initialise it - if it's larger than the old
      ! memory block.

      SELECT CASE (rstorageNode%idataType)
      CASE (ST_SINGLE)
        ALLOCATE(rstorageNode%p_Fsingle1D(isize))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_SINGLE2BYTES,DP)
      CASE (ST_DOUBLE)
        ALLOCATE(rstorageNode%p_Ddouble1D(isize))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_DOUBLE2BYTES,DP)
      CASE (ST_INT)
        ALLOCATE(rstorageNode%p_Iinteger1D(isize))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_INT2BYTES,DP)
      CASE (ST_LOGICAL)
        ALLOCATE(rstorageNode%p_Blogical1D(isize))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_LOGICAL2BYTES,DP)
      CASE (ST_CHAR)
        ALLOCATE(rstorageNode%p_Schar1D(isize))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_CHAR2BYTES,DP)
      CASE DEFAULT
        PRINT *,'Error: unknown mem type'
        CALL sys_halt()
      END SELECT

      IF (isize > isizeOld) &
        CALL storage_initialiseNode (rstorageNode,cinitNewBlock,isizeOld+1_I32)

      ! Copy old data?
      IF (bcopyData) THEN
        isizeCopy=MIN(isize,isizeOld)
        SELECT CASE (rstorageNode%idataType)
        CASE (ST_SINGLE)
          CALL lalg_copyVectorSngl (p_rnode%p_Fsingle1D(1:isizeCopy),&
                                    rstorageNode%p_Fsingle1D(1:isizeCopy))
        CASE (ST_DOUBLE)
          CALL lalg_copyVectorDble (p_rnode%p_Ddouble1D(1:isizeCopy),&
                                    rstorageNode%p_Ddouble1D(1:isizeCopy))
        CASE (ST_INT)
          CALL lalg_copyVectorInt (p_rnode%p_Iinteger1D(1:isizeCopy),&
                                  rstorageNode%p_Iinteger1D(1:isizeCopy))
        CASE (ST_LOGICAL)
          ! Copy by hand
          DO i = 1, isizeCopy
            rstorageNode%p_Blogical1D(i) = p_rnode%p_Blogical1D(i)
          END DO
        CASE (ST_CHAR)
          ! Copy by hand
          DO i = 1, isizeCopy
            rstorageNode%p_Schar1D(i) = p_rnode%p_Schar1D(i)
          END DO
        END SELECT
      END IF

    CASE (2)

      ! Get the size of the old storage node.
      SELECT CASE (p_rnode%idataType)
      CASE (ST_SINGLE)
        Isize2Dold = SHAPE(p_rnode%p_Fsingle2D)
      CASE (ST_DOUBLE)
        Isize2Dold = SHAPE(p_rnode%p_Ddouble2D)
      CASE (ST_INT)
        Isize2Dold = SHAPE(p_rnode%p_Iinteger2D)
      CASE (ST_LOGICAL)
        Isize2Dold = SHAPE(p_rnode%p_Blogical2D)
      CASE (ST_CHAR)
        Isize2Dold = SHAPE(p_rnode%p_Schar2D)
      END SELECT

      ! Do we really have to change anything?
      IF (isize == Isize2Dold(2)) RETURN

      ! Allocate new memory and initialise it - if it's larger than the old
      ! memory block.
      SELECT CASE (rstorageNode%idataType)
      CASE (ST_SINGLE)
        ALLOCATE(rstorageNode%p_Fsingle2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_SINGLE2BYTES,DP)
      CASE (ST_DOUBLE)
        ALLOCATE(rstorageNode%p_Ddouble2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_DOUBLE2BYTES,DP)
      CASE (ST_INT)
        ALLOCATE(rstorageNode%p_Iinteger2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_INT2BYTES,DP)
      CASE (ST_LOGICAL)
        ALLOCATE(rstorageNode%p_Blogical2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_LOGICAL2BYTES,DP)
      CASE (ST_CHAR)
        ALLOCATE(rstorageNode%p_Schar2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_CHAR2BYTES,DP)
      CASE DEFAULT
        PRINT *,'Error: unknown mem type'
        CALL sys_halt()
      END SELECT

      IF (isize > Isize2Dold(2)) &
        CALL storage_initialiseNode (rstorageNode,cinitNewBlock,&
                                     Isize2Dold(2)+1_I32)

      ! Copy old data?
      IF (bcopyData) THEN

        ! Here it's easier than in storage_copy as we can be sure, source and
        ! destination array have the same type!
        SELECT CASE (rstorageNode%idataType)
        CASE (ST_DOUBLE)
          ! Copy by hand
          DO j=1,MIN(SIZE(rstorageNode%p_Ddouble2D,2),Isize2DOld(2))
            DO i=1,SIZE(rstorageNode%p_Ddouble2D,1)
              rstorageNode%p_Ddouble2D(i,j) = p_rnode%p_Ddouble2D(i,j)
            END DO
          END DO

        CASE (ST_SINGLE)
          ! Copy by hand
          DO j=1,MIN(SIZE(rstorageNode%p_Fsingle2D,2),Isize2DOld(2))
            DO i=1,SIZE(rstorageNode%p_Fsingle2D,1)
              rstorageNode%p_Fsingle2D(i,j) = p_rnode%p_Fsingle2D(i,j)
            END DO
          END DO

        CASE (ST_INT)
          ! Copy by hand
          DO j=1,MIN(SIZE(rstorageNode%p_Iinteger2D,2),Isize2DOld(2))
            DO i=1,SIZE(rstorageNode%p_Iinteger2D,1)
              rstorageNode%p_Iinteger2D(i,j) = p_rnode%p_Iinteger2D(i,j)
            END DO
          END DO

        CASE (ST_LOGICAL)
          ! Copy by hand
          DO j=1,MIN(SIZE(rstorageNode%p_Blogical2D,2),Isize2DOld(2))
            DO i=1,SIZE(rstorageNode%p_Blogical2D,1)
              rstorageNode%p_Blogical2D(i,j) = p_rnode%p_Blogical2D(i,j)
            END DO
          END DO

        CASE (ST_CHAR)
          ! Copy by hand
          DO j=1,MIN(SIZE(rstorageNode%p_Schar2D,2),Isize2DOld(2))
            DO i=1,SIZE(rstorageNode%p_Schar2D,1)
              rstorageNode%p_Schar2D(i,j) = p_rnode%p_Schar2D(i,j)
            END DO
          END DO
        END SELECT

      END IF

    CASE DEFAULT
      PRINT *, 'Error in storage_realloc: Handle ',ihandle, &
               ' is neither 1- nor 2- dimensional!'
      CALL sys_halt()

    END SELECT

    ! Respect also the temporary memory in the total amount of memory used.
    IF ((p_rheap%dtotalMem + rstorageNode%dmemBytes) .GT. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem + rstorageNode%dmemBytes

    ! Release old data
    IF (ASSOCIATED(p_rnode%p_Fsingle1D))  DEALLOCATE(p_rnode%p_Fsingle1D)
    IF (ASSOCIATED(p_rnode%p_Ddouble1D))  DEALLOCATE(p_rnode%p_Ddouble1D)
    IF (ASSOCIATED(p_rnode%p_Iinteger1D)) DEALLOCATE(p_rnode%p_Iinteger1D)
    IF (ASSOCIATED(p_rnode%p_Blogical1D)) DEALLOCATE(p_rnode%p_Blogical1D)
    IF (ASSOCIATED(p_rnode%p_Schar1D))    DEALLOCATE(p_rnode%p_Schar1D)
    IF (ASSOCIATED(p_rnode%p_Fsingle2D))  DEALLOCATE(p_rnode%p_Fsingle2D)
    IF (ASSOCIATED(p_rnode%p_Ddouble2D))  DEALLOCATE(p_rnode%p_Ddouble2D)
    IF (ASSOCIATED(p_rnode%p_Iinteger2D)) DEALLOCATE(p_rnode%p_Iinteger2D)
    IF (ASSOCIATED(p_rnode%p_Blogical2D)) DEALLOCATE(p_rnode%p_Blogical2D)
    IF (ASSOCIATED(p_rnode%p_Schar2D))    DEALLOCATE(p_rnode%p_Schar2D)

    ! Correct the memory statistics
    p_rheap%dtotalMem = p_rheap%dtotalMem &
                      - p_rnode%dmemBytes + rstorageNode%dmemBytes
    IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem

    ! Replace the old node by the new one, finish
    p_rnode = rstorageNode

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_reallocFixed (scall, ilbound, iubound, ihandle, &
             cinitNewBlock, bcopy, rheap)

!<description>
  ! This routine reallocates an existing memory block wih a new desired
  ! size. In case of a multiple-dimension block, the last dimension
  ! is changed. isize is the size of the new memory block / the new size
  ! of the last dimension.
  !
  ! Warning: Reallocation of an array destroys all pointers associated with
  ! the corresponding handle!
!</description>

!<input>

  ! Name of the calling routine
  CHARACTER(LEN=*), INTENT(IN) :: scall

  ! Requested lower bound for the memory block / the new lower bound of the last
  ! dimension in the memory block identified by ihandle
  INTEGER(I32), INTENT(IN) :: ilbound

  ! Requested upper bound for the memory block / the new upper bound of the last
  ! dimension in the memory block identified by ihandle
  INTEGER(I32), INTENT(IN) :: iubound

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO,
  ! ST_NEWBLOCK_NOINIT, ST_NEWBLOCK_ORDERED).
  ! Specifies how to initialise memory if isize > original array size.
  INTEGER, INTENT(IN) :: cinitNewBlock

  ! OPTIONAL: Copy old data.
  ! =TRUE: Copy data of old array to the new one.
  ! =FALSE: Reallocate memory, don't copy old data.
  ! If not specified, TRUE is assumed.
  LOGICAL, INTENT(IN), OPTIONAL :: bcopy

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

  ! Handle of the memory block.
  INTEGER, INTENT(INOUT) :: ihandle

!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    TYPE(t_storageBlock), POINTER :: p_rheap
    TYPE(t_storageNode), POINTER :: p_rnode

    ! New storage node
    TYPE(t_storageNode) :: rstorageNode

    ! size of the new 1-dimensional array
    INTEGER(I32) :: isize

    ! size of the old 1-dimensional array
    INTEGER(I32) :: isizeOld

    ! lower bound of the old 1-dimensional array
    INTEGER(I32) :: ilboundOld

    ! upper bound of the old 1-dimensional array
    INTEGER(I32) :: iuboundOld

    ! lower bound of the 1-dimensional array to be copied
    INTEGER :: ilboundCopy

    ! upper bound of the 1-dimensional array to be copied
    INTEGER :: iuboundCopy

    ! size of the old 2-dimensional array
    INTEGER(I32), DIMENSION(2) :: Isize2Dold

    ! lower bound of the old 2-dimensional array
    INTEGER(I32), DIMENSION(2) :: ilbound2Dold

    ! upper bound of the old 2-dimensional array
    INTEGER(I32), DIMENSION(2) :: iubound2Dold

    INTEGER(I32) :: i,j

    LOGICAL :: bcopyData

    isize=iubound-ilbound+1
    IF (isize .EQ. 0) THEN
      ! Ok, not much to do...
      CALL storage_free(ihandle,rheap)
      RETURN
    END IF

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbase
    END IF

    ! Copy old data?

    IF (PRESENT(bcopy)) THEN
      bcopyData = bcopy
    ELSE
      bcopyData = .TRUE.
    END IF

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Copy the data of the old storage node to rstorageNode. That way
    ! we prepare a new storage node and will replace the old.
    rstorageNode = p_rnode

    ! Are we 1D or 2D?
    SELECT CASE(p_rnode%idimension)

    CASE (1)

      ! Get the size and bounds of the old storage node.
      SELECT CASE (p_rnode%idataType)
      CASE (ST_SINGLE)
        isizeOld   = SIZE(p_rnode%p_Fsingle1D)
        ilboundOld = LBOUND(p_rnode%p_Fsingle1D,1)
        iuboundOld = UBOUND(p_rnode%p_Fsingle1D,1)
      CASE (ST_DOUBLE)
        isizeOld   = SIZE(p_rnode%p_Ddouble1D)
        ilboundOld = LBOUND(p_rnode%p_Ddouble1D,1)
        iuboundOld = UBOUND(p_rnode%p_Ddouble1D,1)
      CASE (ST_INT)
        isizeOld   = SIZE(p_rnode%p_Iinteger1D)
        ilboundOld = LBOUND(p_rnode%p_Iinteger1D,1)
        iuboundOld = UBOUND(p_rnode%p_Iinteger1D,1)
      CASE (ST_LOGICAL)
        isizeOld   = SIZE(p_rnode%p_Blogical1D)
        ilboundOld = LBOUND(p_rnode%p_Blogical1D,1)
        iuboundOld = UBOUND(p_rnode%p_Blogical1D,1)
      CASE (ST_CHAR)
        isizeOld   = SIZE(p_rnode%p_Schar1D)
        ilboundOld = LBOUND(p_rnode%p_Schar1D,1)
        iuboundOld = UBOUND(p_rnode%p_Schar1D,1)
      END SELECT

      ! Do we really have to change anything?
      IF ((ilbound == ilboundOld) .AND. &
          (iubound == iuboundOld)) RETURN

      ! Allocate new memory and initialise it - if it's larger than the old
      ! memory block.

      SELECT CASE (rstorageNode%idataType)
      CASE (ST_SINGLE)
        ALLOCATE(rstorageNode%p_Fsingle1D(ilbound:iubound))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_SINGLE2BYTES,DP)
      CASE (ST_DOUBLE)
        ALLOCATE(rstorageNode%p_Ddouble1D(ilbound:iubound))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_DOUBLE2BYTES,DP)
      CASE (ST_INT)
        ALLOCATE(rstorageNode%p_Iinteger1D(ilbound:iubound))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_INT2BYTES,DP)
      CASE (ST_LOGICAL)
        ALLOCATE(rstorageNode%p_Blogical1D(ilbound:iubound))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_LOGICAL2BYTES,DP)
      CASE (ST_CHAR)
        ALLOCATE(rstorageNode%p_Schar1D(ilbound:iubound))
        rstorageNode%dmemBytes = REAL(isize,DP)*REAL(ST_CHAR2BYTES,DP)
      CASE DEFAULT
        PRINT *,'Error: unknown mem type'
        CALL sys_halt()
      END SELECT

      IF (iubound > iuboundOld) &
          CALL storage_initialiseNode (rstorageNode,cinitNewBlock,iuboundOld+1_I32)
      IF (ilbound < ilboundOld) &
          CALL storage_initialiseNode (rstorageNode,cinitNewBlock,ilbound,ilboundOld-1_I32)

      ! Copy old data?
      IF (bcopyData) THEN
        ilboundCopy=MAX(ilbound,ilboundOld)
        iuboundCopy=MIN(iubound,iuboundOld)
        SELECT CASE (rstorageNode%idataType)
        CASE (ST_SINGLE)
          CALL lalg_copyVectorSngl (p_rnode%p_Fsingle1D(ilboundCopy:iuboundCopy),&
                                    rstorageNode%p_Fsingle1D(ilboundCopy:iuboundCopy))
        CASE (ST_DOUBLE)
          CALL lalg_copyVectorDble (p_rnode%p_Ddouble1D(ilboundCopy:iuboundCopy),&
                                    rstorageNode%p_Ddouble1D(ilboundCopy:iuboundCopy))
        CASE (ST_INT)
          CALL lalg_copyVectorInt (p_rnode%p_Iinteger1D(ilboundCopy:iuboundCopy),&
                                  rstorageNode%p_Iinteger1D(ilboundCopy:iuboundCopy))
        CASE (ST_LOGICAL)
          ! Copy by hand
          DO i = ilboundCopy, iuboundCopy
            rstorageNode%p_Blogical1D(i) = p_rnode%p_Blogical1D(i)
          END DO
        CASE (ST_CHAR)
          ! Copy by hand
          DO i = ilboundCopy, iuboundCopy
            rstorageNode%p_Schar1D(i) = p_rnode%p_Schar1D(i)
          END DO
        END SELECT
      END IF

    CASE (2)

      ! Get the size of the old storage node.
      SELECT CASE (p_rnode%idataType)
      CASE (ST_SINGLE)
        Isize2Dold = SHAPE(p_rnode%p_Fsingle2D)
        ilbound2Dold = LBOUND(p_rnode%p_Fsingle2D)
        iubound2Dold = UBOUND(p_rnode%p_Fsingle2D)
      CASE (ST_DOUBLE)
        Isize2Dold = SHAPE(p_rnode%p_Ddouble2D)
        ilbound2Dold = LBOUND(p_rnode%p_Ddouble2D)
        iubound2Dold = UBOUND(p_rnode%p_Ddouble2D)
      CASE (ST_INT)
        Isize2Dold = SHAPE(p_rnode%p_Iinteger2D)
        ilbound2Dold = LBOUND(p_rnode%p_Iinteger2D)
        iubound2Dold = UBOUND(p_rnode%p_Iinteger2D)
      CASE (ST_LOGICAL)
        Isize2Dold = SHAPE(p_rnode%p_Blogical2D)
        ilbound2Dold = LBOUND(p_rnode%p_Blogical2D)
        iubound2Dold = UBOUND(p_rnode%p_Blogical2D)
      CASE (ST_CHAR)
        Isize2Dold = SHAPE(p_rnode%p_Schar2D)
        ilbound2Dold = LBOUND(p_rnode%p_Schar2D)
        iubound2Dold = UBOUND(p_rnode%p_Schar2D)
      END SELECT

      ! Do we really have to change anything?
      IF ((ilbound == ilbound2Dold(2)) .AND. &
          (iubound == iubound2Dold(2))) RETURN

      ! Allocate new memory and initialise it - if it's larger than the old
      ! memory block.
      SELECT CASE (rstorageNode%idataType)
      CASE (ST_SINGLE)
        ALLOCATE(rstorageNode%p_Fsingle2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            &ilbound:iubound))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_SINGLE2BYTES,DP)
      CASE (ST_DOUBLE)
        ALLOCATE(rstorageNode%p_Ddouble2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            ilbound:iubound))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_DOUBLE2BYTES,DP)
      CASE (ST_INT)
        ALLOCATE(rstorageNode%p_Iinteger2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            ilbound:iubound))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_INT2BYTES,DP)
      CASE (ST_LOGICAL)
        ALLOCATE(rstorageNode%p_Blogical2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            ilbound:iubound))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_LOGICAL2BYTES,DP)
      CASE (ST_CHAR)
        ALLOCATE(rstorageNode%p_Schar2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            ilbound:iubound))
        rstorageNode%dmemBytes = &
             REAL(Isize2Dold(1),DP)*REAL(isize,DP)*REAL(ST_CHAR2BYTES,DP)
      CASE DEFAULT
        PRINT *,'Error: unknown mem type'
        CALL sys_halt()
      END SELECT

      IF (iubound > Iubound2Dold(2)) &
        CALL storage_initialiseNode (rstorageNode,cinitNewBlock,&
                                     Iubound2Dold(2)+1_I32)
      IF (ilbound < Ilbound2Dold(2)) &
        CALL storage_initialiseNode (rstorageNode,cinitNewBlock,&
                                     ilbound,Ilbound2Dold(2)-1_I32)

      ! Copy old data?
      IF (bcopyData) THEN

        ! Here it's easier than in storage_copy as we can be sure, source and
        ! destination array have the same type!
        SELECT CASE (rstorageNode%idataType)
        CASE (ST_DOUBLE)
          ! Copy by hand
          DO j=MAX(ilbound,Ilbound2Dold(2)),MIN(iubound,Iubound2Dold(2))
            DO i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Ddouble2D(i,j) = p_rnode%p_Ddouble2D(i,j)
            END DO
          END DO

        CASE (ST_SINGLE)
          ! Copy by hand
          DO j=MAX(ilbound,Ilbound2Dold(2)),MIN(iubound,Iubound2Dold(2))
            DO i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Fsingle2D(i,j) = p_rnode%p_Fsingle2D(i,j)
            END DO
          END DO

        CASE (ST_INT)
          ! Copy by hand
          DO j=MAX(ilbound,Ilbound2Dold(2)),MIN(iubound,Iubound2Dold(2))
            DO i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Iinteger2D(i,j) = p_rnode%p_Iinteger2D(i,j)
            END DO
          END DO

        CASE (ST_LOGICAL)
          ! Copy by hand
          DO j=MAX(ilbound,Ilbound2Dold(2)),MIN(iubound,Iubound2Dold(2))
            DO i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Blogical2D(i,j) = p_rnode%p_Blogical2D(i,j)
            END DO
          END DO

        CASE (ST_CHAR)
          ! Copy by hand
          DO j=MAX(ilbound,Ilbound2Dold(2)),MIN(iubound,Iubound2Dold(2))
            DO i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Schar2D(i,j) = p_rnode%p_Schar2D(i,j)
            END DO
          END DO
        END SELECT

      END IF

    CASE DEFAULT
      PRINT *, 'Error in storage_realloc: Handle ',ihandle, &
               ' is neither 1- nor 2- dimensional!'
      CALL sys_halt()

    END SELECT

    ! Respect also the temporary memory in the total amount of memory used.
    IF ((p_rheap%dtotalMem + rstorageNode%dmemBytes) .GT. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem + rstorageNode%dmemBytes

    ! Release old data
    IF (ASSOCIATED(p_rnode%p_Fsingle1D))  DEALLOCATE(p_rnode%p_Fsingle1D)
    IF (ASSOCIATED(p_rnode%p_Ddouble1D))  DEALLOCATE(p_rnode%p_Ddouble1D)
    IF (ASSOCIATED(p_rnode%p_Iinteger1D)) DEALLOCATE(p_rnode%p_Iinteger1D)
    IF (ASSOCIATED(p_rnode%p_Blogical1D)) DEALLOCATE(p_rnode%p_Blogical1D)
    IF (ASSOCIATED(p_rnode%p_Schar1D))    DEALLOCATE(p_rnode%p_Schar1D)
    IF (ASSOCIATED(p_rnode%p_Fsingle2D))  DEALLOCATE(p_rnode%p_Fsingle2D)
    IF (ASSOCIATED(p_rnode%p_Ddouble2D))  DEALLOCATE(p_rnode%p_Ddouble2D)
    IF (ASSOCIATED(p_rnode%p_Iinteger2D)) DEALLOCATE(p_rnode%p_Iinteger2D)
    IF (ASSOCIATED(p_rnode%p_Blogical2D)) DEALLOCATE(p_rnode%p_Blogical2D)
    IF (ASSOCIATED(p_rnode%p_Schar2D))    DEALLOCATE(p_rnode%p_Schar2D)

    ! Correct the memory statistics
    p_rheap%dtotalMem = p_rheap%dtotalMem &
                      - p_rnode%dmemBytes + rstorageNode%dmemBytes
    IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem

    ! Replace the old node by the new one, finish
    p_rnode = rstorageNode

  END SUBROUTINE

!************************************************************************

!<function>

  FUNCTION storage_isEqual (ihandle1, ihandle2, rheap1, rheap2) RESULT (bisequal)

!<description>

  ! This routine checks if the content of two different handles is equal

!</description>

!<input>

  ! The first handle
  INTEGER, INTENT(IN) :: ihandle1

  ! The second handle
  INTEGER, INTENT(IN) :: ihandle2

  ! OPTIONAL: first local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap1

  ! OPTIONAL: second local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap2

!</input>

!<result>

  ! .TRUE. if the content is equal.
  LOGICAL :: bisequal

!</result>

!</function>

  ! local variables

  ! Pointer to the heaps
  TYPE(t_storageBlock), POINTER :: p_rheap1,p_rheap2

  ! Identifier for data type
  INTEGER :: idatatype1, idatatype2

  ! Identifier for data dimension
  INTEGER :: idimension1, idimension2

  ! Identifier for size
  INTEGER(I32) :: isize1, isize2
  INTEGER(I32) :: Isize2D1(2), Isize2D2(2)

  ! Auxiliary arrays
  REAL(DP), DIMENSION(:,:), POINTER     :: p_Ddouble2D1,p_Ddouble2D2
  REAL(DP), DIMENSION(:),   POINTER     :: p_Ddouble1D1,p_Ddouble1D2
  REAL(SP), DIMENSION(:,:), POINTER     :: p_Fsingle2D1,p_Fsingle2D2
  REAL(SP), DIMENSION(:),   POINTER     :: p_Fsingle1D1,p_Fsingle1D2
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_Iinteger2D1,p_Iinteger2D2
  INTEGER(I32), DIMENSION(:),   POINTER :: p_Iinteger1D1,p_Iinteger1D2
  LOGICAL, DIMENSION(:,:), POINTER      :: p_Blogical2D1,p_Blogical2D2
  LOGICAL, DIMENSION(:),   POINTER      :: p_Blogical1D1,p_Blogical1D2
  CHARACTER, DIMENSION(:,:), POINTER    :: p_Schar2D1,p_Schar2D2
  CHARACTER, DIMENSION(:),   POINTER    :: p_Schar1D1,p_Schar1D2
  
  ! Get the heaps to use - local or global one.

  IF(PRESENT(rheap1)) THEN
    p_rheap1 => rheap1
  ELSE
    p_rheap1 => rbase
  END IF

  IF(PRESENT(rheap2)) THEN
    p_rheap2 => rheap2
  ELSE
    p_rheap2 => rbase
  END IF

  ! Determine data type
  CALL storage_getdatatype(ihandle1, idatatype1, p_rheap1)
  CALL storage_getdatatype(ihandle2, idatatype2, p_rheap2)

  IF (idatatype1 .NE. idatatype2) THEN
    bisequal = .FALSE.
    RETURN
  END IF

  ! Determine dimension
  CALL storage_getdimension(ihandle1, idimension1, p_rheap1)
  CALL storage_getdimension(ihandle2, idimension2, p_rheap2)

  IF (idimension1 .NE. idimension2) THEN
    bisequal = .FALSE.
    RETURN
  END IF

  ! What dimension do we have?
  SELECT CASE(idimension1)
  CASE (1)

    ! Determine size
    CALL storage_getsize1d(ihandle1, isize1, p_rheap1)
    CALL storage_getsize1d(ihandle2, isize2, p_rheap2)

    IF (isize1 .NE. isize2) THEN
      bisequal = .FALSE.
      RETURN
    END IF

    ! What data type do we have?
    SELECT CASE(idatatype1)
    CASE (ST_INT)
      CALL storage_getbase_int(ihandle1, p_Iinteger1D1, p_rheap1)
      CALL storage_getbase_int(ihandle2, p_Iinteger1D2, p_rheap2)

      bisequal = ALL(p_Iinteger1D1 .EQ. p_Iinteger1D2)

    CASE (ST_SINGLE)
      CALL storage_getbase_single(ihandle1, p_Fsingle1D1, p_rheap1)
      CALL storage_getbase_single(ihandle2, p_Fsingle1D2, p_rheap2)
      
      bisequal = ALL(p_Fsingle1D1 .EQ. p_Fsingle1D2)

    CASE (ST_DOUBLE)
      CALL storage_getbase_double(ihandle1, p_Ddouble1D1, p_rheap1)
      CALL storage_getbase_double(ihandle2, p_Ddouble1D2, p_rheap2)

      bisequal = ALL(p_Ddouble1D1 .EQ. p_Ddouble1D2)
      
    CASE (ST_LOGICAL)
      CALL storage_getbase_logical(ihandle1, p_Blogical1D1, p_rheap1)
      CALL storage_getbase_logical(ihandle2, p_Blogical1D2, p_rheap2)

      bisequal = ALL(p_Blogical1D1 .AND. p_Blogical1D2)
      
    CASE (ST_CHAR)
      CALL storage_getbase_char(ihandle1, p_Schar1D1, p_rheap1)
      CALL storage_getbase_char(ihandle2, p_Schar1D2, p_rheap2)

      bisequal = ALL(p_Schar1D1 .EQ. p_Schar1D2)
    END SELECT

  CASE (2)

    ! Determine size
    CALL storage_getsize2d(ihandle1, Isize2D1, p_rheap1)
    CALL storage_getsize2d(ihandle2, Isize2D2, p_rheap2)

    IF (ANY(Isize2d1 .NE. Isize2D2)) THEN
      bisequal = .FALSE.
      RETURN
    END IF

    ! What data type do we have?
    SELECT CASE(idatatype1)
    CASE (ST_INT)
      CALL storage_getbase_int2d(ihandle1, p_Iinteger2D1, p_rheap1)
      CALL storage_getbase_int2d(ihandle2, p_Iinteger2D2, p_rheap2)

      bisequal = ALL(p_Iinteger2D1 .EQ. p_Iinteger2D2)

    CASE (ST_SINGLE)
      CALL storage_getbase_single2d(ihandle1, p_Fsingle2D1, p_rheap1)
      CALL storage_getbase_single2d(ihandle2, p_Fsingle2D2, p_rheap2)
      
      bisequal = ALL(p_Fsingle2D1 .EQ. p_Fsingle2D2)

    CASE (ST_DOUBLE)
      CALL storage_getbase_double2d(ihandle1, p_Ddouble2D1, p_rheap1)
      CALL storage_getbase_double2d(ihandle2, p_Ddouble2D2, p_rheap2)

      bisequal = ALL(p_Ddouble2D1 .EQ. p_Ddouble2D2)
      
    CASE (ST_LOGICAL)
      CALL storage_getbase_logical2d(ihandle1, p_Blogical2D1, p_rheap1)
      CALL storage_getbase_logical2d(ihandle2, p_Blogical2D2, p_rheap2)

      bisequal = ALL(p_Blogical2D1 .AND. p_Blogical2D2)
      
    CASE (ST_CHAR)
      CALL storage_getbase_char2d(ihandle1, p_Schar2D1, p_rheap1)
      CALL storage_getbase_char2d(ihandle2, p_Schar2D2, p_rheap2)

      bisequal = ALL(p_Schar2D1 .EQ. p_Schar2D2)
    END SELECT

  CASE DEFAULT
    PRINT *, 'Error in storage_realloc: Handle ',ihandle1, &
        ' is neither 1- nor 2- dimensional!'
    CALL sys_halt()
  END SELECT
  END FUNCTION storage_isEqual

END MODULE
