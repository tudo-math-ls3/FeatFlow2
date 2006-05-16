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
!# </purpose>
!##############################################################################

MODULE storage

  USE fsystem
  
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
    REAL(DP) :: dmemBytes = 0_DP
    
    ! Pointer to 1D real array or NULL() if not assigned
    REAL(SP), DIMENSION(:), POINTER       :: p_Ssingle1D   => NULL()

    ! Pointer to 1D double precision array or NULL() if not assigned
    REAL(DP), DIMENSION(:), POINTER       :: p_Ddouble1D   => NULL()

    ! Pointer to 1D integer array or NULL() if not assigned
    INTEGER(I32), DIMENSION(:), POINTER   :: p_Iinteger1D   => NULL()

    ! Pointer to 2D real array or NULL() if not assigned
    REAL(SP), DIMENSION(:,:), POINTER     :: p_Ssingle2D   => NULL()

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
    REAL(DP) :: dtotalMem = 0_DP

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

  !</inputouotput>
  
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
    p_Rdescriptors = rheap%p_Rdescriptors
    p_IfreeHandles = rheap%p_IfreeHandles
    
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
  INTEGER, INTENT(IN) :: ihandle
  
  !</input>

  !<inputoutput>
  
  ! The heap structure where to release the handle from.
  TYPE(t_storageBlock), INTENT(INOUT) :: rheap
  
  !</inputoutput>
  
!</subroutine>

  TYPE(t_storageNode), POINTER :: p_rnode
  
  ! Where is the descriptor of the handle?
  p_rnode => rheap%p_Rdescriptors(ihandle)
  
  ! Clear the descriptor structure
  p_rnode%idataType = ST_NOHANDLE
  p_rnode%idimension = 0
  p_rnode%dmemBytes = 0_DP
  NULLIFY(p_rnode%p_Ssingle1D)
  NULLIFY(p_rnode%p_Ddouble1D)
  NULLIFY(p_rnode%p_Iinteger1D)
  NULLIFY(p_rnode%p_Ssingle2D)
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

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  INTEGER, INTENT(IN) :: cinitNewBlock

  !</input>
  
  !<inputoutput>
  
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

  !</inputouotput>

  !<output>

  ! Handle of the memory block.
  INTEGER :: ihandle

  !</output>
  
!</subroutine>
  
  ! Pointer to the heap 
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode
  
  IF (isize .EQ. 0) THEN
    PRINT *,'storage_new1D Warning: isize=0'
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
  p_rnode%idimension = 1
  p_rnode%sname = sname
  
  ! Allocate memory according to isize:
  
  SELECT CASE (ctype)
  CASE (ST_SINGLE)
    ALLOCATE(p_rnode%p_Ssingle1D(isize))
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
  
  ! Clear the vector if necessary
  IF (cinitNewBlock .EQ. ST_NEWBLOCK_ZERO) THEN
    SELECT CASE (ctype)
    CASE (ST_SINGLE)
      p_rnode%p_Ssingle1D = 0.0_SP
    CASE (ST_DOUBLE)
      p_rnode%p_Ddouble1D = 0.0_DP
    CASE (ST_INT)
      p_rnode%p_Iinteger1D = 0_I32
    END SELECT
  END IF

  END SUBROUTINE


!************************************************************************

!<subroutine>

  SUBROUTINE storage_new2D (scall, sname, isize, ctype, ihandle, &
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
  INTEGER, DIMENSION(2), INTENT(IN) :: isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT)
  INTEGER, INTENT(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  INTEGER, INTENT(IN) :: cinitNewBlock

  !</input>
  
  !<inputoutput>
  
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

  !</inputouotput>

  !<output>

  ! Handle of the memory block.
  INTEGER :: ihandle

  !</output>
  
!</subroutine>
  
  ! Pointer to the heap 
  TYPE(t_storageBlock), POINTER :: p_rheap
  TYPE(t_storageNode), POINTER :: p_rnode
  
  IF ((isize(1) .EQ. 0) .OR. (isize(2) .EQ. 0)) THEN
    PRINT *,'storage_new2D Warning: isize=0'
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
  
  ! Allocate memory according to isize:
  
  SELECT CASE (ctype)
  CASE (ST_SINGLE)
    ALLOCATE(p_rnode%p_Ssingle2D(isize(1),isize(2)))
    p_rnode%dmemBytes = p_rnode%dmemBytes + REAL(isize(1),DP)*REAL(isize(1),DP)*REAL(ST_SINGLE2BYTES)
  CASE (ST_DOUBLE)
    ALLOCATE(p_rnode%p_Ddouble2D(isize(1),isize(2)))
    p_rnode%dmemBytes = p_rnode%dmemBytes + REAL(isize(1),DP)*REAL(isize(1),DP)*REAL(ST_DOUBLE2BYTES)
  CASE (ST_INT)
    ALLOCATE(p_rnode%p_Iinteger2D(isize(1),isize(2)))
    p_rnode%dmemBytes = p_rnode%dmemBytes + REAL(isize(1),DP)*REAL(isize(1),DP)*REAL(ST_INT2BYTES)
  CASE DEFAULT
    PRINT *,'Error: unknown mem type'
    STOP
  END SELECT
  
  p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
  
  ! Clear the vector if necessary
  IF (cinitNewBlock .EQ. ST_NEWBLOCK_ZERO) THEN
    SELECT CASE (ctype)
    CASE (ST_SINGLE)
      p_rnode%p_Ssingle2D = 0.0_SP
    CASE (ST_DOUBLE)
      p_rnode%p_Ddouble2D = 0.0_DP
    CASE (ST_INT)
      p_rnode%p_Iinteger2D = 0_I32
    END SELECT
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

  !</inputouotput>

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
  
  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  ! Release the memory assigned to that handle.
  IF (ASSOCIATED(p_rnode%p_Ssingle1D))  DEALLOCATE(p_rnode%p_Ssingle1D)
  IF (ASSOCIATED(p_rnode%p_Ddouble1D))  DEALLOCATE(p_rnode%p_Ddouble1D)
  IF (ASSOCIATED(p_rnode%p_Iinteger1D)) DEALLOCATE(p_rnode%p_Iinteger1D)
  IF (ASSOCIATED(p_rnode%p_Ssingle2D))  DEALLOCATE(p_rnode%p_Ssingle2D)
  IF (ASSOCIATED(p_rnode%p_Ddouble2D))  DEALLOCATE(p_rnode%p_Ddouble2D)
  IF (ASSOCIATED(p_rnode%p_Iinteger2D)) DEALLOCATE(p_rnode%p_Iinteger2D)
  
  ! Release the handle itself.
  CALL storage_releasehandle (ihandle,p_rheap)

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
    PRINT *,'Wrong handle'
    STOP
  END IF
  
  IF (p_rheap%p_Rdescriptors(ihandle)%idataType .NE. ST_SINGLE) THEN
    PRINT *,'storage_getbase_single: Wrong data format!'
    STOP
  END IF

  ! Get the pointer  
  
  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Ssingle1D
  
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
    PRINT *,'Wrong handle'
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
    PRINT *,'Wrong handle'
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
    PRINT *,'Wrong handle'
    STOP
  END IF
  
  ! Get the pointer  
  
  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Ssingle2D
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE storage_getbase_double2D (ihandle, p_Darray, rheap)
  
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
    PRINT *,'Wrong handle'
    STOP
  END IF
  
  ! Get the pointer  
  
  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D
  
  END SUBROUTINE
  
END MODULE
