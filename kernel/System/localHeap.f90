!CRAP!!!

!##############################################################################
!# ****************************************************************************
!# <name> localHeap </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a small memory-management component.
!#
!# When solving problems on different levels, it's usually necessary to
!# allocate and deallocate memory very often, especially on lower levels.
!# These memory blocks used as temporary memory are usually quite small.
!# This leads to signigicant performance lost due to the frequent
!# allocations and deallocations.
!#
!# This module offers a possibility to encounter such performance loss
!# by allocating one large array of memory, so called 'local heap'. 
!# The application can request memory of this block for local usage and 
!# release it when it's no more needed.
!#
!# The basic memory management supports the following operations:
!# 1.) Create temporary memory
!#     Creates a temporary memory block of specific size with a specific
!#     number of 'handles', i.e. identifiers for the subblocks.
!# 2.) Request memory
!#     This returns a handle and a pointer to a memory block in the
!#     global array
!# 3.) Release memory
!#     Releases the memory block specified by the handle
!# 4.) Release temporary memory
!# In the basic memory management, it's adviseable to release the memory
!# blocks in EXACTLY THE OPPOSITE ORDER as they were allocated - 
!# otherwise holes will be created in these 'local heaps', leading
!# to performance loss.
!#
!# If more memory is requested than available, the memory is allocated
!# on the global heap. If no memory is available there, an error
!# will be returned.
!# If more handles are requested than available, the handle list will
!# be extended automatically.
!#
!# The heap supports only Integer and Real data blocks at the moment.
!# </purpose>
!##############################################################################

MODULE localHeap

  USE fsystem
  
  IMPLICIT NONE

!<constants>

  !<constantblock description="Kind values for data that can be stored.">

  ! kind value for integer data
  INTEGER, PARAMETER :: LOCHEAP_PRECINT  = I32

  ! kind value for standard integers
  INTEGER, PARAMETER :: LOCHEAP_PRECREAL = DP
  
  ! kind value for indices used for indexing memory blocks. I32 allowes
  ! to access memory blocks with ~2.000.000.000 entries, I64 more.
  INTEGER, PARAMETER :: LOCHEAP_PRECIDX  = I32
  
  !</constantblock>

  !<constantblock description="General constants for temporary memory management.">

  ! Maximum size of arrays in the temporary memory management.
  ! If a the memory is exceeded, the standard ALLOCATE is used.
  INTEGER(LOCHEAP_PRECIDX), PARAMETER :: LOCHEAP_MAXTEMPMEM  = HUGE(1_LOCHEAP_PRECIDX)/4_LOCHEAP_PRECIDX

  !</constantblock>
  
!</constants>

!<types>

  !<typeblock>

  ! PRIVATE STRUCTURE: configures how to map a handle to a memory block
  ! in a data array.

  TYPE t_localHeapNode
    PRIVATE 
      
    ! Type of memory associated to a handle. 0=real, 1=integer
    INTEGER(INT) :: imemType = 0
    
    ! Starting address of the memory block in an array.
    ! =0, if the handle to this node is not in use
    ! =-1 if the memory is allocated on the heap with ALLOCATE
    INTEGER(LOCHEAP_PRECIDX) :: istartIndex = 0
    
    ! Ending address of the memory block in an array
    ! =0, if the handle to this node is not in use
    ! =-1 if the memory is allocated on the heap with ALLOCATE
    INTEGER(LOCHEAP_PRECIDX) :: iendIndex   = 0
    
    ! Pointer to real data array if this node refers to real data.
    ! NULL=not assigned.
    INTEGER(LOCHEAP_PRECREAL), DIMENSION(:), POINTER :: p_Ddata => NULL()

    ! Pointer to real data array if this node refers to integer data
    ! NULL=not assigned.
    INTEGER(LOCHEAP_PRECINT), DIMENSION(:), POINTER :: p_Idata => NULL()
    
  END TYPE

  !</typeblock>
  
  !<typeblock>
  
  ! The basic structure for temporary memory management.
  
  TYPE t_localHeap
  
    ! Size of memory in use in integer-array
    INTEGER(LOCHEAP_PRECIDX)                         :: nimemInUse = 0
    
    ! A pointer to an allocated memory block for integer data -
    ! or NULL if not allocated.
    INTEGER(LOCHEAP_PRECINT), DIMENSION(:), POINTER :: p_Idata  => NULL()
    
    ! Size of memory in use in real-array
    INTEGER(LOCHEAP_PRECIDX)                         :: ndmemInUse = 0
    
    ! A pointer to an allocated memory block for real data -
    ! or NULL if not allocated.
    REAL(LOCHEAP_PRECREAL), DIMENSION(:), POINTER    :: p_Ddata => NULL()
    
    ! Handle-array. When k a handle, the IstartAddress(k) configures
    ! whether the handle is associated to a real-subarray in Ddata,
    ! or an integer-array in Idata. It furthermore configures the
    ! starting and ending address of the array in Idata/Ddata.
    TYPE(t_localHeapNode), DIMENSION(:), POINTER :: p_RhandleDef => NULL()
    
    ! Number of handles available in p_RhandleDef
    INTEGER(INT)                                  :: nhandlesAvailable = 0
    
    ! Number of handles that are in use; <= SIZE(p_RhandleDef).
    INTEGER(INT)                                  :: nhandlesInUse = 0
    
    ! Number of holes if handles < maximum handle is released
    INTEGER(INT)                                  :: nholes = 0

  END TYPE
  
  !</typeblock>

!</types>

  INTERFACE locheap_localHeapAllocate
    MODULE PROCEDURE locheap_localHeapAllocateReal
    MODULE PROCEDURE locheap_localHeapAllocateInt
  END INTERFACE

  INTERFACE locheap_localHeapRelease
    MODULE PROCEDURE locheap_localHeapReleaseReal
    MODULE PROCEDURE locheap_localHeapReleaseInt
  END INTERFACE

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE locheap_localHeapInitialise (localHeap, ninteger, nreal, &
                                          nhandles)
                                       
  !<description>
  
  ! This routine initialises a local-heap structure. It allocates memory
  ! for nintegers integer-values, nreal double values and nhandles handles
  ! on the global heap.
  
  !</description>
  
  !<input>
  
  ! Number of integer values to be reserved on the heap
  INTEGER(LOCHEAP_PRECIDX), INTENT(IN) :: ninteger

  ! Number of real values to be reserved on the heap
  INTEGER(LOCHEAP_PRECIDX), INTENT(IN) :: nreal
  
  ! Number of handles to be reserved on the heap
  INTEGER(INT), INTENT(IN)          :: nhandles
  
  !</inout>
  
  !<output>
  
  ! The local-heap structure to be initialised
  TYPE(t_localHeap), INTENT(OUT) :: localHeap
  
  !</output>

!</subroutine>

  IF ( (nhandles .EQ. 0) .OR. ((ninteger .EQ. 0) .AND. (nreal .EQ. 0)) ) THEN
    PRINT *,'Error: No memory block to allocate?!?'
    STOP
  END IF

  ! Reserve memory:
  
  ALLOCATE(localHeap%p_RhandleDef(nhandles))
  localHeap%nhandlesAvailable = nhandles
  localHeap%nhandlesInUse     = 0
  
  IF (ninteger .NE. 0) THEN
    ALLOCATE(localHeap%p_Idata(ninteger))
  ELSE
    NULLIFY(localHeap%p_Idata)
  END IF

  IF (nreal .NE. 0) THEN
    ALLOCATE(localHeap%p_Ddata(nreal))
  ELSE
    NULLIFY(localHeap%p_Ddata)
  END IF
  
  ! That's it.

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE locheap_localHeapAllocate (localHeap, nsize, ihandle)
  
  !<description>
  
  ! This routine tries to allocate isize elements on the 'local' heap
  ! identified by localHeap. If this does not work, the memory is
  ! allocated on the global heap.
  
  !</description>
  
  !<inputoutput>

  ! The t_localHeap structure identifying the local heap
  TYPE(t_localHeap), INTENT(INOUT) :: localHeap

  !</inputoutput>

  !<input>
  
  ! The length of the array to be allocated
  INTEGER(LOCHEAP_PRECIDX), INTENT(IN)    :: nsize
  
  !</input>
  
  !<output>
  
  ! The handle of the memory block that is allocated.
  INTEGER(INT), INTENT(OUT)            :: ihandle
  
  !</output>
  
!</subroutine>

  ! local variables
  INTEGER(LOCHEAP_PRECIDX)    :: nsizeleft, istartidx, iendidx
  INTEGER(INT)                :: i,j,iisLocal

  IF (.NOT. ASSOCIATED(localHeap%p_RhandleDef)) THEN
    PRINT *,'Error: local heap not initialised'
    STOP
  END IF

  IF (.NOT. ASSOCIATED(localHeap%p_Idata)) THEN
    PRINT *,'Error: local heap not initialised for integers'
    STOP
  END IF
  
  ! No more handles available?
  IF (localHeap%nhandlesInUse .GE. localHeap%nhandlesAvailable) THEN
  
    iisLocal = NO
  
    ! Do we have holes in the heap? If yes, we might be lucky to find
    ! a hole large enough to hold our data!
    IF (localHeap%nholes .GT. 0) THEN
      
      ! Try to find the first hole that is large enough to hold the memory
      ! block. As nhandlesInUse >= nhandlesAvailable, the last handle
      ! is allocated for sure!
      istartidx = 1
      iendidx = 1
      i = 1
      DO WHILE (i <= localHeap%nhandlesAvailable-1)
        
        IF (localHeap%p_RhandleDef(i)%istartIndex .EQ. 0) THEN
        
          ! Find the length of the hole - at least the last node exists!
          ! (so the loop sets iendidx for sure)
          DO j = i+1,localHeap%nhandlesAvailable
            IF (localHeap%p_RhandleDef(j)%istartIndex .NE. 0) THEN
              iendidx = localHeap%p_RhandleDef(j)%istartIndex - 1
              iisLocal = YES
              EXIT
            END IF
          END DO
          
          IF (iendidx-istartidx .GT. nsize) THEN
            ! Yeah, we can use the hole!
            localHeap%p_RhandleDef(ihandle)%istartIndex = istartidx
            localHeap%p_RhandleDef(ihandle)%iendIndex   = iendidx
            localHeap%p_RhandleDef(ihandle)%imemType    = 1
            iisLocal = YES
            EXIT
          END IF
          
          ! Hole not large enough - continue the search at the next handle
          ! that is not a hole!
          i = j
          
        ELSE
          
          i = i+1
          
        END IF
      
      END DO
      
      ! No hole available - use the global heap!
    
    END IF
    
    IF (iisLocal .EQ. NO) THEN
    ! Allocate on the global heap
    localHeap%nhandlesInUseOutOfBounds = localHeap%nhandlesInUseOutOfBounds+1
    ihandle = localHeap%nhandlesInUse
    ALLOCATE(p_Iarray(nsize))
    RETURN
  
  END IF
  
            p_Iarray => localHeap%p_Idata( &
                            localHeap%p_RhandleDef(ihandle)%istartIndex : &
                            localHeap%p_RhandleDef(ihandle)%iendIndex)




  ! Memory full?
  
  nsizeleft = SIZE(localHeap%p_Idata) - localHeap%nimemInUse 
  
  IF ( (nsize .GT. SIZE(localHeap%p_Idata)) .OR. &
       (nsize .GT. nsizeleft) ) THEN

    ! Allocate on the global heap
    localHeap%nhandlesInUse = localHeap%nhandlesInUse+1
    ihandle = localHeap%nhandlesInUse
    ALLOCATE(p_Iarray(nsize))
    
    ! Mark in the handle structure that the memory was allocated globally...
    
    localHeap%p_RhandleDef(ihandle)%istartIndex = -1
    localHeap%p_RhandleDef(ihandle)%iendIndex   = -1
    localHeap%p_RhandleDef(ihandle)%imemType    = 1
    RETURN
    
  END IF
  
  ! Ok, we can use our memory management.
  ! Allocate the memory in the 'local' heap.

  localHeap%nhandlesInUse = localHeap%nhandlesInUse+1
  ihandle = localHeap%nhandlesInUse
  localHeap%p_RhandleDef(ihandle)%istartIndex = localHeap%nimemInUse+1
  localHeap%p_RhandleDef(ihandle)%iendIndex   = localHeap%nimemInUse+nsize
  localHeap%p_RhandleDef(ihandle)%imemType    = 1
  localHeap%nimemInUse = localHeap%nimemInUse+nsize
  
  p_Iarray => localHeap%p_Idata( &
                  localHeap%p_RhandleDef(ihandle)%istartIndex : &
                  localHeap%p_RhandleDef(ihandle)%iendIndex)
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE locheap_localHeapAllocateReal (localHeap, nsize, ihandle, &
                                            p_Darray)
  
  !<description>
  
  ! This routine tries to allocate isize reals on the 'local' heap
  ! identified by localHeap. If this does not work, the memory is
  ! allocated on the global heap.
  
  !</description>
  
  !<inputoutput>

  ! The t_localHeap structure identifying the local heap
  TYPE(t_localHeap), INTENT(INOUT) :: localHeap

  !</inputoutput>

  !<input>
  
  ! The length of the array to be allocated
  INTEGER(LOCHEAP_PRECIDX), INTENT(IN)    :: nsize
  
  !</input>
  
  !<output>
  
  ! The handle of the memory block that is allocated.
  INTEGER(INT), INTENT(OUT)            :: ihandle
  
  ! A pointer to the memory block
  REAL(LOCHEAP_PRECREAL),DIMENSION(:),POINTER,INTENT(OUT) :: p_Darray
  
  !</output>
  
!</subroutine>
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE locheap_localHeapReleaseInt (localHeap, ihandle, p_Iarray)
  
  !<description>
  
  ! This routine releases the memory associated to a handle from the
  ! 'local' heap.
  ! Remark: Memory blocks must be released in exactly the
  ! opposite order as they were allocated!
  
  !</description>
  
  !<inputoutput>

  ! The t_localHeap structure identifying the local heap
  TYPE(t_localHeap), INTENT(INOUT) :: localHeap

  ! The handle of the memory block that should be released.
  ! Will be set to 0.
  INTEGER(INT), INTENT(INOUT)             :: ihandle

  ! A pointer to the memory block associated to ihandle.
  ! Is nullified.
  INTEGER(LOCHEAP_PRECINT),DIMENSION(:),POINTER,INTENT(INOUT) :: p_Iarray

  !</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(INT) :: i

  IF (.NOT. ASSOCIATED(localHeap%p_RhandleDef)) THEN
    PRINT *,'Error: local heap not initialised'
    STOP
  END IF

  IF (.NOT. ASSOCIATED(localHeap%p_Idata)) THEN
    PRINT *,'Error: local heap not initialised for integers'
    STOP
  END IF
  
  ! Handle valid?
  IF ( (ihandle .LT. 0) .OR. &
       ( (ihandle .GT. localHeap%nhandlesInUse) .AND. &
         (localHeap%nhandlesInUse < localHeap%nhandlesAvailable) ) ) THEN
    PRINT *,'Error: Invalid handle'
    STOP
  END IF
  
  ! If the handle is 'out of bounds' concerning the handle-array, release it from
  ! the global heap.
  IF ( ihandle .GT. localHeap%nhandlesAvailable ) THEN
    
    DEALLOCATE(p_Iarray)
    
    ! Decrease the number of global arrays in use out of the bounds
    ! of the handle-array
    
    localHeap%nhandlesInUseOutOfBounds = localHeap%nhandlesInUseOutOfBounds-1
    ihandle = 0
    RETURN
    
  END IF

  ! Is the handle the last one?
  IF ( (ihandle .EQ. localHeap%nhandlesInUse) ) THEN
  
    ! Easier, quicker case, release the memory - either from the local or
    ! from the global heap.
  
    ! If the memory was allocated on the global heap, release ir
    IF ( localHeap%p_RhandleDef(ihandle)%istartIndex .EQ. -1 ) THEN
      DEALLOCATE(p_Iarray)
      localHeap%nhandlesInUse = localHeap%nhandlesInUse-1
      ihandle = 0
      RETURN
    END IF
    
    ! If it was allocated on the local heap, release it from the array.
    NULLIFY(p_Iarray)
    localHeap%nhandlesInUse = localHeap%nhandlesInUse-1
    localHeap%p_RhandleDef(ihandle)%istartIndex = 0
    localHeap%p_RhandleDef(ihandle)%iendIndex = 0
    ihandle = 0

    ! Now we have to figure out how many memory is still allocated on the 
    ! 'local' heap. This removes 'holes' which might appear by releasing
    ! handles inbetween.
    ! Go back from the last known handle till we find one with iendIndex<>0.
    ! That's the last position occupied in the local heap.
    
    DO i=localHeap%nhandlesInUse,1,-1
      IF (localHeap%p_RhandleDef(ihandle)%iendIndex.NE.0) EXIT
    END DO
    
    IF (I.NE.0) THEN
      
      ! Ok, now we know the amount of available memory
      localHeap%nimemInUse = localHeap%p_RhandleDef(ihandle)%iendIndex
    
    ELSE
    
      ! Memory is completely free
      localHeap%nimemInUse = 0
    
    END IF

  ELSE
  
    ! The handle is not the last one. Holes might come up!
  
    ! If the memory was allocated on the global heap, release ir
    IF ( localHeap%p_RhandleDef(ihandle)%istartIndex .EQ. -1 ) THEN
      DEALLOCATE(p_Iarray)
      ! Don't decrease the number of handles in use as we released
      ! memory 'in the inner'. nhandlesInUse is always the maximum
      ! number of handles in use, not counting holes.
      localHeap%p_RhandleDef(ihandle)%istartIndex = 0
      localHeap%p_RhandleDef(ihandle)%iendIndex = 0
      ihandle = 0
      RETURN
    END IF
    
    ! Release the memory from the 'local' heap
    ! If it was allocated on the local heap, release it from the array.
    NULLIFY(p_Iarray)
    ihandle = 0

    ! Don't decrease the number of handles in use as we released
    ! memory 'in the inner'. nhandlesInUse is always the maximum
    ! number of handles in use, not counting holes.
    !
    ! Don't change the amount of available memory, as this creates a
    ! 'hole' in the memory array.

  END IF
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE locheap_localHeapReleaseReal (localHeap, ihandle, p_Darray)
  
  !<description>
  
  ! This routine releases the memory associated to a handle from the
  ! 'local' heap.
  ! Remark: Memory blocks must be released in exactly the
  ! opposite order as they were allocated!
  
  !</description>
  
  !<inputoutput>

  ! The t_localHeap structure identifying the local heap
  TYPE(t_localHeap), INTENT(INOUT) :: localHeap

  ! The handle of the memory block that should be released.
  ! Will be set to 0.
  INTEGER(INT), INTENT(INOUT)             :: ihandle

  ! A pointer to the memory block associated to ihandle.
  ! Is nullified.
  REAL(LOCHEAP_PRECREAL),DIMENSION(:),POINTER,INTENT(INOUT) :: p_Darray

  !</inputoutput>

!</subroutine>

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE locheap_localHeapDone(localHeap,ifreeOnPurpose)
                                       
  !<description>
  
  ! This routine cleans up a local-heap structure. All allocated memory
  ! is released from the global heap - except for those memory blocks
  ! allocated on the global heap!
    
  !</description>
  
  !<inputoutput>
  
  ! The local-heap structure to be initialised
  TYPE(t_localHeap), INTENT(INOUT) :: localHeap
  
  !</inputoutput>

  !<input>

  ! Optional. Skip memory leak check. If set to YES, there will be no
  ! warning message if the heap is released and there are still
  ! allocated memory blocks in the heap. All memory belonging to
  ! the local heap (also the memory allocated on the global heap if
  ! necessary) will be released without warnings.
  INTEGER(INT), INTENT(IN)         :: ifreeOnPurpose
  
  !</input>
  
!</subroutine>

  ! local variables
  INTEGER(INT) :: inoWarn, i, iWarn

  ! Release memory

  IF (.NOT. ASSOCIATED(localHeap%p_RhandleDef)) THEN
    PRINT *,'Error: local heap not initialised'
    STOP
  END IF
  
  inoWarn = NO
  IF (PRESENT(ifreeOnPurpose)) inoWarn = ifreeOnPurpose
  
  IF (localHeap%nhandlesInUse .NE. 0) THEN
  
    IF ((inoWarn .EQ. NO) .AND. ) THEN
      PRINT *,'Warning: Memory leak! Releasing forgotten memory...'
    END IF
    
    ! Loop through the heap, release all handles that are still associated
    
    iWarn = 0
    DO i=1,localHeap%
    
    
    
  END IF
  
  
  IF (ASSOCIATED(localHeap%p_Ddata)) THEN
    DEALLOCATE(localHeap%p_Ddata)
  END IF
  
  IF (ASSOCIATED(localHeap%p_Idata)) THEN
    DEALLOCATE(localHeap%p_Idata)
  END IF

  IF (ASSOCIATED(localHeap%p_RhandleDef)) THEN
    DEALLOCATE(localHeap%p_RhandleDef)
  END IF

  ! Clean up, finish

  localHeap%nimemInUse = 0
  localHeap%ndmemInUse = 0
  localHeap%nhandlesAvailable = 0
  localHeap%nhandlesInUse = 0

  END SUBROUTINE

  ! ***************************************************************************

END MODULE
