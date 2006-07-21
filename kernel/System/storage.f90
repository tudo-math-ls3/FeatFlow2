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
!# PRECISION floating point variables as well as I32 integer variables.
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
!#      -> Determine pointer associated to a handle for singles, doubles
!#         or 32-Bit integers
!#
!#  7.) storage_getbase_single2D,
!#      storage_getbase_double2D,
!#      storage_getbase_int2D,
!#      -> Determine pointer associated to a handle for singles, doubles
!#         or 32-Bit integers, 2D array
!#
!#  8.) storage_copy
!#      -> Copies the content of one array to another.
!#
!#  9.) storage_clear
!#      -> Clears an array by overwriting the entries with 0.
!#
!# 10.) storage_getsize = storage_getsize1d / storage_getsize2d
!#      -> Get the length of an array on the heap.
!# 
!# </purpose>
!##############################################################################

MODULE storage

  USE fsystem
  USE linearalgebra
  
  IMPLICIT NONE
  
!<constants>

  !<constantblock description="Storage block type identifiers">
  
  !defines an non-allocated storage block handle
  INTEGER, PARAMETER :: ST_NOHANDLE = 0

  !storage block contains single floats
  INTEGER, PARAMETER :: ST_SINGLE = 1

  !storage block contains double floats
  INTEGER, PARAMETER :: ST_DOUBLE = 2

  !storage block contains ints
  INTEGER, PARAMETER :: ST_INT = 3
  
  !</constantblock>

  !<constantblock description="Constants for initialisation of memory on allocation">
  
  ! init new storage block with zeros
  INTEGER, PARAMETER :: ST_NEWBLOCK_ZERO = 0

  ! no init new storage block
  INTEGER, PARAMETER :: ST_NEWBLOCK_NOINIT = 1

  ! init new storage block with 1,2,3,...,n
  INTEGER, PARAMETER :: ST_NEWBLOCK_ORDERED = 2
  
  !</constantblock>
  
  !<constantblock description="Constants for calculating memory">

  ! How many bytes has an integer?
  INTEGER :: ST_INT2BYTES = I32

  ! How many bytes has a standard real?
  INTEGER :: ST_SINGLE2BYTES = SP

  ! How many bytes has a double precision real?
  INTEGER :: ST_DOUBLE2BYTES = DP

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
    REAL(SP), DIMENSION(:), POINTER       :: p_Fsingle1D   => NULL()

    ! Pointer to 1D double precision array or NULL() if not assigned
    REAL(DP), DIMENSION(:), POINTER       :: p_Ddouble1D   => NULL()

    ! Pointer to 1D integer array or NULL() if not assigned
    INTEGER(I32), DIMENSION(:), POINTER   :: p_Iinteger1D   => NULL()

    ! Pointer to 2D real array or NULL() if not assigned
    REAL(SP), DIMENSION(:,:), POINTER     :: p_Fsingle2D   => NULL()

    ! Pointer to 2D double precision array or NULL() if not assigned
    REAL(DP), DIMENSION(:,:), POINTER     :: p_Ddouble2D   => NULL()

    ! Pointer to 2D integer array or NULL() if not assigned
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_Iinteger2D  => NULL()

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
    
    ! Maximum amount of memory that was in use ofer the whole lifetime
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
    MODULE PROCEDURE storage_new2D
  END INTERFACE

  INTERFACE storage_getsize
    MODULE PROCEDURE storage_getsize1D
    MODULE PROCEDURE storage_getsize2D
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
    STOP
  END IF
  
  ! Handles available?
  
  IF (rheap%ihandlesInUse .GE. rheap%nhandlesTotal) THEN
  
    ! All handles are in use. We have to modify our ring to accept more
    ! andles.
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
  NULLIFY(p_rnode%p_Fsingle2D)
  NULLIFY(p_rnode%p_Ddouble2D)
  NULLIFY(p_rnode%p_Iinteger2D)
  
  ! Handle ihandle is available now - put it to the list of available handles.
  rheap%p_ilastFreeHandle = MOD(rheap%p_ilastFreeHandle,rheap%nhandlesTotal) + 1
  rheap%p_IfreeHandles (rheap%p_ilastFreeHandle) = ihandle
   
  rheap%ihandlesInUse = rheap%ihandlesInUse - 1
  
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
  INTEGER, INTENT(IN) :: isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT)
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
  INTEGER :: ihandle

!</output>
  
!</subroutine>
  
  ! Pointer to the heap 
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode

  ! variable for ordering 1,2,3,...,N
  INTEGER :: iorder
  
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
  
  ! Get a new handle
  ihandle = storage_newhandle (p_rheap)

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Initialise the content
  
  p_rnode%idataType = ctype
  p_rnode%idimension = 1
  p_rnode%sname = sname
  
  ! Allocate memory according to isize:
  
  SELECT CASE (ctype)
  CASE (ST_SINGLE)
    ALLOCATE(p_rnode%p_Fsingle1D(isize))
    p_rnode%dmemBytes = p_rnode%dmemBytes + REAL(isize,DP)*REAL(ST_SINGLE2BYTES)
  CASE (ST_DOUBLE)
    ALLOCATE(p_rnode%p_Ddouble1D(isize))
    p_rnode%dmemBytes = p_rnode%dmemBytes + REAL(isize,DP)*REAL(ST_DOUBLE2BYTES)
  CASE (ST_INT)
    ALLOCATE(p_rnode%p_Iinteger1D(isize))
    p_rnode%dmemBytes = p_rnode%dmemBytes + REAL(isize,DP)*REAL(ST_INT2BYTES)
  CASE DEFAULT
    PRINT *,'Error: unknown mem type'
    STOP
  END SELECT
  
  p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
  IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
    p_rheap%dtotalMemMax = p_rheap%dtotalMem
  
  ! Clear the vector if necessary
  IF (cinitNewBlock .EQ. ST_NEWBLOCK_ZERO) THEN
    SELECT CASE (ctype)
    CASE (ST_SINGLE)
      p_rnode%p_Fsingle1D = 0.0_SP
    CASE (ST_DOUBLE)
      p_rnode%p_Ddouble1D = 0.0_DP
    CASE (ST_INT)
      p_rnode%p_Iinteger1D = 0_I32
    END SELECT
  END IF

  ! Impose ordering 1,2,3,...,N if necessary
  IF (cinitNewBlock .EQ. ST_NEWBLOCK_ORDERED) THEN
    SELECT CASE (ctype)
    CASE (ST_SINGLE)
      DO iorder=1,isize
        p_rnode%p_Fsingle1D(iorder) = REAL(iorder,SP)
      END DO
    CASE (ST_DOUBLE)
      DO iorder=1,isize
        p_rnode%p_Ddouble1D(iorder) = REAL(iorder,DP)
      END DO
    CASE (ST_INT)
      DO iorder=1,isize
        p_rnode%p_Iinteger1D(iorder) = INT(iorder,I32)
      END DO
    END SELECT
  END IF

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
  INTEGER, DIMENSION(2), INTENT(IN) :: Isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT)
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
  INTEGER :: ihandle

!</output>
  
!</subroutine>
  
  ! Pointer to the heap 
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode
  
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
  
  ! Get a new handle
  ihandle = storage_newhandle (p_rheap)

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Initialise the content
  
  p_rnode%idataType = ctype
  p_rnode%idimension = 2
  p_rnode%sname = sname
  
  ! Allocate memory according to Isize:
  
  SELECT CASE (ctype)
  CASE (ST_SINGLE)
    ALLOCATE(p_rnode%p_Fsingle2D(Isize(1),Isize(2)))
    p_rnode%dmemBytes = p_rnode%dmemBytes + &
                        REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_SINGLE2BYTES)
  CASE (ST_DOUBLE)
    ALLOCATE(p_rnode%p_Ddouble2D(Isize(1),Isize(2)))
    p_rnode%dmemBytes = p_rnode%dmemBytes + &
                        REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_DOUBLE2BYTES)
  CASE (ST_INT)
    ALLOCATE(p_rnode%p_Iinteger2D(Isize(1),Isize(2)))
    p_rnode%dmemBytes = p_rnode%dmemBytes + &
                        REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_INT2BYTES)
  CASE DEFAULT
    PRINT *,'Error: unknown mem type'
    STOP
  END SELECT
  
  p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
  IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
    p_rheap%dtotalMemMax = p_rheap%dtotalMem
  
  ! Clear the vector if necessary
  IF (cinitNewBlock .EQ. ST_NEWBLOCK_ZERO) THEN
    SELECT CASE (ctype)
    CASE (ST_SINGLE)
      p_rnode%p_Fsingle2D = 0.0_SP
    CASE (ST_DOUBLE)
      p_rnode%p_Ddouble2D = 0.0_DP
    CASE (ST_INT)
      p_rnode%p_Iinteger2D = 0_I32
    END SELECT
  END IF

  IF (cinitNewBlock .EQ. ST_NEWBLOCK_ORDERED) THEN
    PRINT *, "Error: ordering not available for multidimensional array"
    STOP
  END IF
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
    STOP
  END IF
  
  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)
  
  ! Is the node associated at all?
  IF (p_rnode%idataType .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_free: Trying to release nonexistent handle!'
    PRINT *,'Handle number: ',ihandle
    STOP
  END IF

  ! Release the memory assigned to that handle.
  IF (ASSOCIATED(p_rnode%p_Fsingle1D))  DEALLOCATE(p_rnode%p_Fsingle1D)
  IF (ASSOCIATED(p_rnode%p_Ddouble1D))  DEALLOCATE(p_rnode%p_Ddouble1D)
  IF (ASSOCIATED(p_rnode%p_Iinteger1D)) DEALLOCATE(p_rnode%p_Iinteger1D)
  IF (ASSOCIATED(p_rnode%p_Fsingle2D))  DEALLOCATE(p_rnode%p_Fsingle2D)
  IF (ASSOCIATED(p_rnode%p_Ddouble2D))  DEALLOCATE(p_rnode%p_Ddouble2D)
  IF (ASSOCIATED(p_rnode%p_Iinteger2D)) DEALLOCATE(p_rnode%p_Iinteger2D)
  
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
    STOP
  END IF
  
  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)
  
  ! Is the node associated at all?
  IF (p_rnode%idataType .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_clear: Trying to release nonexistent handle!'
    PRINT *,'Handle number: ',ihandle
    STOP
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
    END SELECT
  CASE (2)
    SELECT CASE (p_rnode%idataType)
    CASE (ST_SINGLE)
      p_rnode%p_Fsingle2D = 0.0_SP
    CASE (ST_DOUBLE)
      p_rnode%p_Ddouble2D = 0.0_DP
    CASE (ST_INT)
      p_rnode%p_Iinteger2D = 0_I32
    END SELECT
  CASE DEFAULT
    PRINT *,'storage_clear: invalid dimension.'
    STOP
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
    STOP
  END IF
  
  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)
  
  ! Is the node associated at all?
  IF (p_rnode%idataType .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getsize1D: Handle invalid!'
    PRINT *,'Handle number: ',ihandle
    STOP
  END IF

  ! What are we?
  IF (p_rnode%idimension .NE. 1) THEN
    PRINT *,'Error in storage_getsize1D: Handle ',ihandle,' is not 1-dimensional!'
    STOP
  END IF
  
  SELECT CASE (p_rnode%idataType)
  CASE (ST_SINGLE)
    isize = SIZE(p_rnode%p_Fsingle1D)
  CASE (ST_DOUBLE)
    isize = SIZE(p_rnode%p_Ddouble1D)
  CASE (ST_INT)
    isize = SIZE(p_rnode%p_Iinteger1D)
  CASE DEFAULT
    PRINT *,'Error in storage_getsize1D: Invalid data type!' 
    STOP
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
    STOP
  END IF
  
  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)
  
  ! Is the node associated at all?
  IF (p_rnode%idataType .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in storage_getsize2D: Handle invalid!'
    PRINT *,'Handle number: ',ihandle
    STOP
  END IF

  ! What are we?
  IF (p_rnode%idimension .NE. 2) THEN
    PRINT *,'Error in storage_getsize1D: Handle ',ihandle,' is not 2-dimensional!'
    STOP
  END IF
  
  SELECT CASE (p_rnode%idataType)
  CASE (ST_SINGLE)
    Isize = SIZE(p_rnode%p_Fsingle2D)
  CASE (ST_DOUBLE)
    Isize = SIZE(p_rnode%p_Ddouble2D)
  CASE (ST_INT)
    Isize = SIZE(p_rnode%p_Iinteger2D)
  CASE DEFAULT
    PRINT *,'Error in storage_getsize2D: Invalid data type!' 
    STOP
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
    STOP
  END IF
  
  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_INT) THEN
    PRINT *,'storage_getbase_int: Wrong data format!'
    STOP
  END IF
  
  ! Get the pointer  
  
  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_single (ihandle, p_Sarray, rheap)
  
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
    STOP
  END IF
  
  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_SINGLE) THEN
    PRINT *,'storage_getbase_single: Wrong data format!'
    STOP
  END IF

  ! Get the pointer  
  
  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_double (ihandle, p_Darray, rheap)
  
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
    STOP
  END IF
  
  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_DOUBLE) THEN
    PRINT *,'storage_getbase_double: Wrong data format!'
    STOP
  END IF

  ! Get the pointer  
  
  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D
  
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
    STOP
  END IF
  
  ! Get the pointer  
  
  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D
  
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
    STOP
  END IF
  
  ! Get the pointer  
  
  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_double2D (ihandle, p_Darray, rheap)
  
!<description>
  
  ! This routine returns the pointer to a handle associated to an
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
    STOP
  END IF
  
  ! Get the pointer  
  
  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D
  
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
      STOP
    END IF
    IF (.NOT. ASSOCIATED(p_rheap%p_Rdescriptors)) THEN
      PRINT *,'storage_copy: Heap not initialised!'
      STOP
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
          CALL storage_new ('storage_copy',p_rsource%sname,SIZE(p_rsource%p_Ddouble1D),&
                            ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_SINGLE)
          CALL storage_new ('storage_copy',p_rsource%sname,SIZE(p_rsource%p_Fsingle1D),&
                            ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        CASE (ST_INT)
          CALL storage_new ('storage_copy',p_rsource%sname,SIZE(p_rsource%p_Iinteger1D),&
                            ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
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
        END SELECT
      END SELECT
      
    END IF
    
    p_rdest => p_rheap%p_Rdescriptors(h_dest)

    ! 1D/2D the same?
    IF (p_rsource%idimension .NE. p_rdest%idimension) THEN
      PRINT *,'storage_copy: Structure different!'
      STOP
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
          STOP
        END SELECT
        
      CASE (ST_SINGLE)
        SELECT CASE (p_rdest%idataType)
        CASE (ST_DOUBLE)
          CALL lalg_copyVectorDblSngl (p_rsource%p_Fsingle1D,p_rdest%p_Ddouble1D)
        CASE (ST_SINGLE)
          CALL lalg_copyVectorSngl (p_rsource%p_Fsingle1D,p_rdest%p_Fsingle1D)
        CASE DEFAULT
          PRINT *,'storage_copy: Unsupported data type combination'
          STOP
        END SELECT
        
      CASE (ST_INT)
        IF (p_rdest%idataType .EQ. ST_INT) THEN
          CALL lalg_copyVectorInt (p_rsource%p_Iinteger1D,p_rdest%p_Iinteger1D)
        ELSE
          PRINT *,'storage_copy: Unsupported data type combination'
          STOP
        END IF
      CASE DEFAULT
        PRINT *,'storage_copy: Unknown data type'
        STOP
      END SELECT  
      
    CASE (2)
      SELECT CASE (p_rsource%idataType)
      CASE (ST_DOUBLE)
        IF ((SIZE(p_rsource%p_Ddouble2D,1) .NE. SIZE(p_rdest%p_Ddouble2D,1)) .OR.&
            (SIZE(p_rsource%p_Ddouble2D,2) .NE. SIZE(p_rdest%p_Ddouble2D,2))) THEN
          PRINT *,'storage_copy: Structure different!'
          STOP
        END IF
        
        SELECT CASE (p_rdest%idataType)
        CASE (ST_DOUBLE)
          ! Copy by hand
          DO j=1,SIZE(p_rsource%p_Ddouble2D,2)
            DO i=1,SIZE(p_rsource%p_Ddouble2D,1)
              p_rdest%p_Ddouble2D(i,j) = p_rsource%p_Ddouble2D(i,j)
            END DO
          END DO
          
        CASE (ST_SINGLE)
          ! Copy by hand
          DO j=1,SIZE(p_rsource%p_Fsingle2D,2)
            DO i=1,SIZE(p_rsource%p_Fsingle2D,1)
              p_rdest%p_Fsingle2D(i,j) = p_rsource%p_Ddouble2D(i,j)
            END DO
          END DO
        
        CASE (ST_INT)
          ! Copy by hand
          DO j=1,SIZE(p_rsource%p_Iinteger2D,2)
            DO i=1,SIZE(p_rsource%p_Iinteger2D,1)
              p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            END DO
          END DO
          
        END SELECT
        
      CASE (ST_SINGLE)
        IF ((SIZE(p_rsource%p_Fsingle2D,1) .NE. SIZE(p_rdest%p_Fsingle2D,1)) .OR.&
            (SIZE(p_rsource%p_Fsingle2D,2) .NE. SIZE(p_rdest%p_Fsingle2D,2))) THEN
          PRINT *,'storage_copy: Structure different!'
          STOP
        END IF

        SELECT CASE (p_rdest%idataType)
        CASE (ST_DOUBLE)
          ! Copy by hand
          DO j=1,SIZE(p_rsource%p_Ddouble2D,2)
            DO i=1,SIZE(p_rsource%p_Ddouble2D,1)
              p_rdest%p_Ddouble2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            END DO
          END DO

        CASE (ST_SINGLE)
          ! Copy by hand
          DO j=1,SIZE(p_rsource%p_Fsingle2D,2)
            DO i=1,SIZE(p_rsource%p_Fsingle2D,1)
              p_rdest%p_Fsingle2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            END DO
          END DO

        CASE (ST_INT)
          ! Copy by hand
          DO j=1,SIZE(p_rsource%p_Iinteger2D,2)
            DO i=1,SIZE(p_rsource%p_Iinteger2D,1)
              p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            END DO
          END DO
        END SELECT

      CASE (ST_INT)
        IF ((SIZE(p_rsource%p_Iinteger2D,1) .NE. SIZE(p_rdest%p_Iinteger2D,1)) .OR.&
            (SIZE(p_rsource%p_Iinteger2D,2) .NE. SIZE(p_rdest%p_Iinteger2D,2))) THEN
          PRINT *,'storage_copy: Structure different!'
          STOP
        END IF
        IF (p_rdest%idataType .NE. ST_INT) THEN
          PRINT *,'storage_copy: unsupported data type combination'
          STOP
        END IF
        
        ! Copy by hand
        DO j=1,SIZE(p_rsource%p_Iinteger2D,2)
          DO i=1,SIZE(p_rsource%p_Iinteger2D,1)
            p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Iinteger2D(i,j)
          END DO
        END DO

      CASE DEFAULT
        PRINT *,'storage_copy: Unknown data type'
        STOP
      END SELECT  
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
    
    PRINT *,'Heap statistics:'
    PRINT *,'----------------'
    IF (PRESENT(bprintHandles)) THEN
      IF (bprintHandles .AND. (p_rheap%ihandlesInUse .GT. 0)) THEN
        PRINT *,'Handles on the heap: '
        PRINT *
        ! Loop through the heap and search allocated handles
        DO i=1,SIZE(p_rheap%p_IfreeHandles)
          IF (p_rheap%p_Rdescriptors(i)%idataType .NE. ST_NOHANDLE) THEN
            IF (p_rheap%p_Rdescriptors(i)%idimension .EQ. 1) THEN
              PRINT *,'Handle ',i,', 1D, Length=',&
                      INT(p_rheap%p_Rdescriptors(i)%dmemBytes,I32),&
                      ', Type=',p_rheap%p_Rdescriptors(i)%idataType,&
                      ' Name=',TRIM(ADJUSTL(p_rheap%p_Rdescriptors(i)%sname))
            ELSE
              PRINT *,'Handle ',i,', 2D, Length=',&
                      INT(p_rheap%p_Rdescriptors(i)%dmemBytes,I32),&
                      ', Type=',p_rheap%p_Rdescriptors(i)%idataType,&
                      ' Name=',TRIM(ADJUSTL(p_rheap%p_Rdescriptors(i)%sname))
            END IF
          END IF
        END DO
        PRINT *
      END IF
    END IF
    
    PRINT *,'Number of allocated handles:     ',p_rheap%ihandlesInUse
    PRINT *,'Current total number of handles: ',SIZE(p_rheap%p_IfreeHandles)
    PRINT *,'Memory in use:                   ',INT(p_rheap%dtotalMem)
    PRINT *,'Maximum used memory:             ',INT(p_rheap%dtotalMemMax)

  END SUBROUTINE
  
END MODULE
