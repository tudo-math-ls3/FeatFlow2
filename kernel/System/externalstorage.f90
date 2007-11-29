!##############################################################################
!# ****************************************************************************
!# <name> externalstorage </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module realises an external storage management for storing memory
!# blocks/arrays on external devices like on a hard disc. The following
!# routines can be found here.
!#
!#  1.) exstor_init
!#      -> Initialises the external storage management
!#
!#  2.) exstor_attachRamdrive
!#      -> Attaches a new RamDrive container to the list of storage containers
!#
!#  3.) exstor_attachDirectory
!#      -> Attaches a file system container to the list of storage containers
!#
!#  4.) exstor_done
!#      -> Cleans up the external storage management
!#
!#  5.) exstor_info
!#      -> Prints statistics about the storage management to the terminal
!#
!#  6.) exstor_new = exstor_new1D/exstor_new2D
!#      -> Creates a new memory block
!#  
!#  7.) exstor_free
!#      -> Releases a memory block 
!#
!#  8.) exstor_getsize = exstor_getsize1D/exstor_getsize2D
!#      -> Returns the size of a memory block
!#
!#  9.) exstor_getdata_int,
!#      exstor_getdata_single,
!#      exstor_getdata_double,
!#      exstor_getdata_storage
!#      -> Reads 1D data from a storage container
!#
!# 10.) exstor_setdata_int
!#      exstor_setdata_single,
!#      exstor_setdata_double,
!#      exstor_setdata_storage
!#      -> Writes 1D data to a storage container
!#
!# 11.) exstor_clear
!#      -> Clears an array by overwriting the entries with 0 (or .FALSE. for LOGICALs).
!#
!# 12.) exstor_copy
!#      -> Copies the content of one array to another.
!#
!# The external storage management consists of a number of 'storage containers'.
!# A storage container can be for example a directory on a hard disc. Each 
!# storage container contains a set of data blocks of different data type.
!#
!# Data blocks are them primarily maintained by the following routines:
!#
!# - exstor_new 
!#   -> will create a new memory block in such a container and identify it
!#      with a handle. 
!#
!# - exstor_getdata_xxxx and exstor_setdata_xxxx
!#   -> allow to read/write from/to a memory block. The block can be released
!#
!# - exstor_free
!#   -> Releases a memory block
!# 
!# It's possible to maintain more than one type of storage container.
!# By default, exstor_init creates at least one storage container, the so
!# called 'RamDrive container' (with Id=1), which is actually a mapping to the
!# standard memory management of the storage.f90. All data stored in this
!# RamDrive container will be stored in the main memory of the computer.
!# Using exstor_attachRamdrive, it's possible to add more RamDrive containers,
!# although this should rarely be necessary.
!#
!# The next type of storage container is the so called 'directory' storage
!# container. Using exstor_attachDirectory, a directory is attached as
!# storage container to the storage management. The data blocks in this
!# container are realised as files on the file system of the computer.
!#
!# By default, exstor_new creates new memory blocks on the last data container
!# attached to the memory management which is large enough to hold the data.
!# However, it's possible to instruct exstor_new to use a special container.
!# The container with id=1 is always the RamDrive container, while the type
!# of the containers with id=2,3,4,... depends the order of creation using
!# the exstor_attachXXXX commands.
!#
!# The typical usage of this module is as follows:
!#
!#   ! Initialise external storage managemeng
!#   CALL exstor_init (999,100)
!#   
!#   ! Append a directory as 2nd storage container
!#   CALL exstor_attachDirectory('./ff2storage')
!#
!#   ! Create data blocks with exstor_new.
!#   ! Read/write data from/to the container using
!#   ! exstor_getdata_XXXX / exstor_setdata_XXXX.
!#   ! Remove the data with exstor_free.
!# 
!#   ! Clean up the external storage management.
!#   CALL exstor_done()
!#
!# </purpose>
!#
!# Missing: Implementation for chars and logicals.
!#
!##############################################################################

MODULE externalstorage

  USE fsystem
  USE genoutput
  USE storage
  USE io

  IMPLICIT NONE

!<constants>

!<constantblock description="Type flags identifying external storage container.">

  ! Undefined container
  INTEGER, PARAMETER :: EXSTOR_CONT_UNDEFINED   = -1

  ! In-memory/RamDrive container. The 'external storage' is actually 'in memory' and
  ! not on any hard disc or similar external storage.
  INTEGER, PARAMETER :: EXSTOR_CONT_RAMDRIVE    = 0
  
  ! The external storage container is a subdirectory on the hard disc.
  INTEGER, PARAMETER :: EXSTOR_CONT_DIRECTORY   = 1

!</constantblock>

!<constantblock description="Strategy constants for choosing a container.">

  ! Always choose the last existing container. If it's full, choose the
  ! last but one, etc.
  INTEGER, PARAMETER :: EXSTOR_STRAT_LASTCONTAINER = 0

!</constantblock>

!</constants>

!<types>

  !<typeblock>

  ! This type block specifies a container for external storage.
  ! This may be for example a directory on a hard disc.

  TYPE t_exstorageContainer

    PRIVATE
    
    ! Type of the external storage container. One of the EXSTOR_CONT_xxxx
    ! flags.
    INTEGER :: ctype = EXSTOR_CONT_UNDEFINED
    
    ! Maximum storage (in Megabytes) that this storage container can handle.
    ! =-1: Infinite memory available.
    INTEGER :: imaxStorageMB = -1
    
    ! Number of bytes currently allocated by this storage container.
    ! This is an integer number encoded as a double precision number
    ! in order to capture even the size of large memory drives.
    REAL(DP) :: dcurrentStorage = 0.0_DP
    
    ! For ctype=EXSTOR_CONT_DIRECTORY: Name of the directory containing
    ! the files with data. The string contains a "/" at the end.
    CHARACTER(LEN=SYS_STRLEN) :: spath = './'
    
    ! For ctype=EXSTOR_CONT_DIRECTORY: Basic filename of the data files.
    CHARACTER(LEN=SYS_NAMELEN) :: sfilename = 'feat2exstorage'
    
    ! For ctype=EXSTOR_CONT_DIRECTORY: Whether to save the data formatted 
    ! or unformatted.
    LOGICAL :: bformatted = .FALSE.

  END TYPE

  !</typeblock>

  !<typeblock>

  ! Type block for describing a handle. This collects the number of the
  ! handle, the storage amount associated to it, the pointer to the memory
  ! location etc.

  TYPE t_exstorageNode
  
    PRIVATE

    ! Type of data associated to the handle (ST_NOHANDLE, ST_SINGLE,
    ! ST_DOUBLE, ST_INT, ST_LOGICAL, ST_CHAR)
    INTEGER :: cdataType = ST_NOHANDLE

    ! Dimension associated to the handle (0=not assigned, 1=1D, 2=2D array)
    INTEGER :: idimension = 0

    ! The name of the array that is associated to that handle
    CHARACTER(LEN=SYS_NAMELEN) :: sname = ''

    ! Amount of memory (in bytes) associated to this block.
    ! We store that as a double to allow storing numbers > 2GB !
    REAL(DP) :: dmemBytes = 0.0_DP
    
    ! Size of the data array. For 1D arrays, the array size can be found
    ! in Isize(1).
    INTEGER(I32), DIMENSION(2) :: Isize = (/0,0/)

    ! Id of the external storage container that contains the data
    INTEGER :: icontainerId = 0
    
    ! Whether the storage block is bound to a container or not.
    ! TRUE=The memory block is bound to container icontainerId.
    ! FALSE=The memory management system may automatically move the memory 
    !       block from one container to another if necessary.
    LOGICAL :: bcontainerBound = .FALSE.

    ! Flag that identifies whether this storage block is initialised
    ! somehow.
    ! 
    ! Flag that specifies whether the data in this vector is 
    ! treated as being completely zero.
    ! ST_NEWBLOCK_NOINIT:  The storage block contains data.
    ! ST_NEWBLOCK_ZERO:    The storage block is to be treated as zero.
    ! ST_NEWBLOCK_ORDERED: The storage block is to be treated as initialised
    !                      by a sequence of numbers (1,2,3,...)
    INTEGER :: cinitNewBlock = ST_NEWBLOCK_NOINIT
    
    ! This is a handle to a memory block. The handle is valid if the
    ! storage block is realised in the main memory as part of a RamDrive
    ! container. (A RamDrive container uses the standard memory management 
    ! of the storage.f90 and saves memory blocks in the main memory.)
    INTEGER :: istorageHandle = ST_NOHANDLE
    
    ! if the storage block is realised as a file in a directory, this
    ! variable contains the filename of the file.
    CHARACTER(LEN=SYS_NAMELEN) :: sfilename = ''
    
  END TYPE

  !</typeblock>

  !<typeblock>

  ! This block represents a heap that maintains singole, double precision
  ! and integer data. It contains a list of t_storageNode elements for all
  ! the handles.
  ! There's one global object of this type for the global storage management,
  ! but if necessary, an algorithm can create such a block locally, too,
  ! to prevent conflicts with the global memory.

  TYPE t_exstorageBlock

    PRIVATE

    ! An array of t_exstorageContainer objects that identifies all
    ! possible storage containers that are handled by this storage
    ! block.
    TYPE(t_exstorageContainer), DIMENSION(:), POINTER :: p_RstorageContainers => NULL()

    ! Strategy how to choose a container in p_RstorageContainers where
    ! to save data. One of the EXSTOR_STRAT_xxxx costants. By default,
    ! this is EXSTOR_STRAT_LASTCONTAINER, so always the last container
    ! providing enough free memory is chosen to save new data.
    INTEGER :: ccontainerStrategy = EXSTOR_STRAT_LASTCONTAINER

    ! An array of t_exstorageNode objects corresponding to the handles.
    ! Each entry identifies a memory block in an external storage container.
    ! Can be dynamically extended if there are not enough handles available.
    TYPE(t_exstorageNode), DIMENSION(:), POINTER :: p_Rdescriptors => NULL()

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
  TYPE(t_exstorageBlock), SAVE, TARGET :: rbaseExternal

!</globals>

  INTERFACE exstor_new
    MODULE PROCEDURE exstor_new1D
    MODULE PROCEDURE exstor_new2D
  END INTERFACE
  
  INTERFACE exstor_getsize
    MODULE PROCEDURE exstor_getsize1D
    MODULE PROCEDURE exstor_getsize2D
  END INTERFACE
  
  PRIVATE :: exstor_newhandle
  PRIVATE :: exstor_releasehandle
  PRIVATE :: exstor_getContainer

CONTAINS

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_init(ihandleCount, ihandlesDelta, rheap)

!<description>
  ! This routine initializes the storage management for external storage.
  ! ihandleCount is the initial number of handles maintained by the
  ! storage routines. If there are not enough free handles, the number
  ! of handles are increased by ihandlesDelta (which is initially set
  ! to 1/2*ihandleCount if not given).
  ! rheap allows to specify a 'local' heap structure to initialise.
  ! If not given, the global external memory management is initialised.
  !
  ! By default, there is exactly one 'RamDrive' container attached to
  ! the storage management, so all memory maintained here is held in the
  ! Ram of the PC. To change this behaviour, attach a hard drive container
  ! to the heap by calling exstor_attachDirectory.
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
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!</subroutine>

    ! local variables

    ! the real 'handle-delta'
    INTEGER :: ihandles, ihDelta

    ! Pointer to the heap to initialise
    TYPE(t_exstorageBlock), POINTER :: p_rheap

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
      p_rheap => rbaseexternal
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
    
    ! Attach a RamDrive container with arbitrary memory size.
    CALL exstor_attachRamdrive(-1,p_rheap)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_attachRamdrive(imaxMem,rheap)

!<description>
  ! This routine attaches a Ramdrive container to the external storage
  ! management.
!</description>

!<input>
  ! OPTIONAL: Maximum size of the RamDrive (in Megabytes).
  ! -1 or not present = arbitrary.
  INTEGER, INTENT(IN), OPTIONAL :: imaxMem
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is initialised.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_exstorageContainer), DIMENSION(:), POINTER :: p_RstorageContainers
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Pointer to the heap to initialise
    TYPE(t_exstorageBlock), POINTER :: p_rheap

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    ! Create a new container.
    IF (.NOT. ASSOCIATED(p_rheap%p_RstorageContainers)) THEN
      ALLOCATE(p_RstorageContainers(1))
    ELSE
      ALLOCATE(p_RstorageContainers(SIZE(p_rheap%p_RstorageContainers)+1))
      p_RstorageContainers(1:SIZE(p_rheap%p_RstorageContainers)) = &
          p_rheap%p_RstorageContainers(:)
      
      DEALLOCATE(p_rheap%p_RstorageContainers)
    END IF
    p_rheap%p_RstorageContainers => p_RstorageContainers
    
    p_rcontainer => p_rheap%p_RstorageContainers(SIZE(p_rheap%p_RstorageContainers))
    
    ! Initialise the container as a RamDrive container
    p_rcontainer%ctype = EXSTOR_CONT_RAMDRIVE
    IF(PRESENT(imaxMem)) p_rcontainer%imaxStorageMB = imaxMem
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_attachDirectory(spath,imaxMem,sfilename,bformatted,rheap)

!<description>
  ! This routine attaches a directory-on-disc container to the external storage
  ! management.
!</description>

!<input>
  ! Path to the directory that saves the data.
  ! The directory must exist. If this is "", the current directory is used.
  CHARACTER(LEN=*), INTENT(IN) :: spath

  ! OPTIONAL: If set to TRUE, the data in the container will be saved
  ! in a formatted file format. If set to FALSE (default), the data will
  ! be saved unformatted (which is machine dependent but faster).
  LOGICAL, INTENT(IN), OPTIONAL :: bformatted
  
  ! OPTIONAL: Maximum size of the container (in Megabytes).
  ! -1 or not present = arbitrary.
  INTEGER, INTENT(IN), OPTIONAL :: imaxMem

  ! OPTIONAL: Basic filename of files that are stored on disc. The files will get
  ! the name "[filename].[handlenr]". If not specified, a default filename
  ! will be used.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sfilename
  
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is initialised.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_exstorageContainer), DIMENSION(:), POINTER :: p_RstorageContainers
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Pointer to the heap to initialise
    TYPE(t_exstorageBlock), POINTER :: p_rheap

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    ! Create a new container.
    IF (.NOT. ASSOCIATED(p_rheap%p_RstorageContainers)) THEN
      ALLOCATE(p_RstorageContainers(1))
    ELSE
      ALLOCATE(p_RstorageContainers(SIZE(p_rheap%p_RstorageContainers)+1))
      p_RstorageContainers(1:SIZE(p_rheap%p_RstorageContainers)) = &
          p_rheap%p_RstorageContainers(:)
      
      DEALLOCATE(p_rheap%p_RstorageContainers)
    END IF
    p_rheap%p_RstorageContainers => p_RstorageContainers
    
    p_rcontainer => p_rheap%p_RstorageContainers(SIZE(p_rheap%p_RstorageContainers))
    
    ! Initialise the container as a Directory container
    p_rcontainer%ctype = EXSTOR_CONT_DIRECTORY
    IF (spath .EQ. '') THEN
      p_rcontainer%spath = './'
    ELSE
      p_rcontainer%spath = TRIM(spath)//'/'
    END IF
    IF(PRESENT(imaxMem)) p_rcontainer%imaxStorageMB = imaxMem
    IF(PRESENT(sfilename)) p_rcontainer%sfilename = sfilename
    IF(PRESENT(bformatted)) p_rcontainer%bformatted = bformatted
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_done(rheap)

!<description>
  ! This routine cleans up the storage management. All data on the
  ! heap is released from memory.
!</description>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is cleaned up.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap to initialise
    TYPE(t_exstorageBlock), POINTER :: p_rheap

    INTEGER :: i,ihandle

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    ! Delete all data from the heap
    DO i = 1,SIZE(p_rheap%p_Rdescriptors)
      ! Don't pass i as handle as storage_free will set the handle
      ! passed to it to 0!
      ihandle = i
      IF (p_rheap%p_Rdescriptors(i)%cdataType .NE. ST_NOHANDLE) &
        CALL exstor_free(ihandle,rheap)
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
    
    ! Release the storage containers
    DEALLOCATE(p_rheap%p_RstorageContainers)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_info(bprintContainers, bprintHandles,rheap)

!<description>
  ! This routine prints information about the current memory consumption
  ! in a memory block to screen.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the information about the storage containers 
  ! is printed to the terminal.
  LOGICAL, INTENT(IN), OPTIONAL :: bprintContainers

  ! OPTIONAL: If set to TRUE, the handles still remaining in the
  ! heap together with their names are printed to the terminal.
  LOGICAL, INTENT(IN), OPTIONAL :: bprintHandles

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_exstorageBlock), INTENT(IN), TARGET, OPTIONAL :: rheap
!</input>

!</subroutine>

  ! local variables
  INTEGER :: i

  ! Pointer to the heap
  TYPE(t_exstorageBlock), POINTER :: p_rheap

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    CALL output_line ('External memory heap statistics:')
    CALL output_line ('--------------------------------')
    IF (PRESENT(bprintHandles)) THEN
      IF (bprintHandles .AND. (p_rheap%ihandlesInUse .GT. 0)) THEN
        CALL output_line ('Handles on the heap: ')
        CALL output_lbrk ()
        ! Loop through the heap and search allocated handles
        DO i=1,SIZE(p_rheap%p_IfreeHandles)
          IF (p_rheap%p_Rdescriptors(i)%cdataType .NE. ST_NOHANDLE) THEN
            IF (p_rheap%p_Rdescriptors(i)%idimension .EQ. 1) THEN
              CALL output_line ( &
                   'Handle ' // TRIM(sys_siL(i,10)) // ', 1D, Length=' // &
                   TRIM(sys_siL(INT(p_rheap%p_Rdescriptors(i)%dmemBytes),15)) //&
                   ', Type=' // TRIM(sys_siL(p_rheap%p_Rdescriptors(i)%cdataType,15)) //&
                   ' Name=' // TRIM(ADJUSTL(p_rheap%p_Rdescriptors(i)%sname)) )
            ELSE
              CALL output_line ( &
                   'Handle ' // TRIM(sys_siL(i,10)) // ', 2D, Length=' // &
                   TRIM(sys_siL(INT(p_rheap%p_Rdescriptors(i)%dmemBytes),15)) // &
                   ', Type=' // TRIM(sys_siL(p_rheap%p_Rdescriptors(i)%cdataType,15)) //&
                   ' Name=' // TRIM(ADJUSTL(p_rheap%p_Rdescriptors(i)%sname)) )
            END IF
          END IF
        END DO
        CALL output_lbrk ()
      END IF
    END IF

    IF (PRESENT(bprintContainers)) THEN
      IF (bprintContainers) THEN
        CALL output_line ('Storage containers: ')
        
        DO i=1,SIZE(p_rheap%p_RstorageContainers)
          CALL output_lbrk ()
          
          CALL output_line ('Storage container:  ' // TRIM(sys_siL(i,10)))
          
          SELECT CASE (p_rheap%p_RstorageContainers(i)%ctype)
          CASE (EXSTOR_CONT_RAMDRIVE)
            CALL output_line ('Type:               RamDrive')
          CASE (EXSTOR_CONT_DIRECTORY)
            CALL output_line ('Type:               Directory ('//&
                TRIM(p_rheap%p_RstorageContainers(i)%spath)//')')
          END SELECT
          
          IF (p_rheap%p_RstorageContainers(i)%imaxStorageMB .EQ. -1) THEN
            CALL output_line ('Max. memory(MB):    infinite')
          ELSE
            CALL output_line ('Max. memory(MB):    '//&
              TRIM(sys_siL(p_rheap%p_RstorageContainers(i)%imaxStorageMB,10)))
          END IF
          
          CALL output_line ('Current memory(MB): '//&
              TRIM(sys_siL(INT(p_rheap%p_RstorageContainers(i)%dcurrentStorage&
                   /1000000.0_DP,I32),10)))

        END DO
      
        CALL output_lbrk ()
      END IF
    END IF
                      
    CALL output_line ('Number of storage containers:    '//&
                      TRIM(sys_siL(SIZE(p_rheap%p_RstorageContainers),15)))
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

!<function>

  INTEGER FUNCTION exstor_newhandle (rheap) RESULT(ihandle)

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
  TYPE(t_exstorageBlock), INTENT(INOUT) :: rheap

!</inputoutput>

!</function>

    ! local variables
    TYPE(t_exstorageNode), DIMENSION(:), POINTER :: p_Rdescriptors => NULL()
    INTEGER, DIMENSION(:), POINTER :: p_IfreeHandles => NULL()
    INTEGER :: i

    IF (rheap%nhandlesTotal .LE. 0) THEN
      CALL output_line ('Heap not initialised!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'exstor_newhandle')
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

  SUBROUTINE exstor_releasehandle (ihandle,rheap)

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
  TYPE(t_exstorageBlock), INTENT(INOUT) :: rheap
!</inputoutput>

!</subroutine>

    TYPE(t_exstorageNode), POINTER :: p_rnode
    INTEGER :: icontainer

    ! Where is the descriptor of the handle?
    p_rnode => rheap%p_Rdescriptors(ihandle)

    ! Subtract the memory amount from the statistics
    rheap%dtotalMem = rheap%dtotalMem - p_rnode%dmemBytes
    
    ! The same for the corresponding container
    icontainer = p_rnode%icontainerId
    rheap%p_RstorageContainers(icontainer)%dcurrentStorage = &
        rheap%p_RstorageContainers(icontainer)%dcurrentStorage - p_rnode%dmemBytes

    ! Clear the descriptor structure
    p_rnode%cdataType = ST_NOHANDLE
    p_rnode%idimension = 0
    p_rnode%dmemBytes = 0.0_DP
    p_rnode%icontainerId = 0
    p_rnode%bcontainerBound = .FALSE.
    p_rnode%cinitNewBlock = ST_NEWBLOCK_NOINIT
    p_rnode%Isize(1:2) = (/0,0/)

    ! Handle ihandle is available now - put it to the list of available handles.
    rheap%p_ilastFreeHandle = MOD(rheap%p_ilastFreeHandle,rheap%nhandlesTotal) + 1
    rheap%p_IfreeHandles (rheap%p_ilastFreeHandle) = ihandle

    rheap%ihandlesInUse = rheap%ihandlesInUse - 1

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_free (ihandle, rheap)

!<description>
  ! This routine releases a handle from a heap and deallocates the
  ! associated memory. ihandle is set to ST_NOHANDLE upon return.
!</description>

!<inputoutput>

  ! Handle of the memory block to be releases
  INTEGER :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!</subroutine>

  ! local variables

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .LE. ST_NOHANDLE) THEN
      CALL output_line ('Releasing ST_NOHANDLE is not allowed!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_free')
      CALL sys_halt()
    END IF

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Is the node associated at all?
    IF (p_rnode%cdataType .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Trying to release nonexistent handle: '//sys_siL(ihandle,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_free')
      CALL sys_halt()
    END IF

    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    ! Release the memory assigned to that handle. How this is released depends
    ! on the container type.
    SELECT CASE(p_rcontainer%ctype)

    CASE (EXSTOR_CONT_RAMDRIVE)
      ! The handle identifies another handle in the main memory storage subsystem.
      ! Note that the handle may be ST_NOHANDLE -- which is the case if no data
      ! is associated to this up to now.
      IF (p_rnode%istorageHandle .NE. ST_NOHANDLE) &
        CALL storage_free (p_rnode%istorageHandle)

    CASE (EXSTOR_CONT_DIRECTORY)
      ! The handle identifies a file in a directory. Delete the file.
      CALL io_deleteFile (TRIM(p_rcontainer%spath)//p_rnode%sfilename)

    END SELECT

    ! Release the handle itself.
    CALL exstor_releasehandle (ihandle,p_rheap)

    ! And finally reset the handle to ST_NOHANDLE.
    ihandle = ST_NOHANDLE

  END SUBROUTINE

!************************************************************************

!<function>

  INTEGER FUNCTION exstor_getContainer (rheap,dsize)

!<description>
  ! This routine determines a container that is large enough to hold
  ! a memory block of size dsize.
!</description>

!<input>
  ! Local heap structure where a container should be searched.
  TYPE(t_exstorageBlock), INTENT(IN), TARGET, OPTIONAL :: rheap
  
  ! Number of bytes (encoded as double) which are to be stored.
  REAL(DP), INTENT(IN) :: dsize
!</input>

!<result>
  ! An id of a container that is large enough to hold the memory.
  ! If no container is found, Id 1 (the standard Id for the RamDrive
  ! container) is returned.
!</result>

!</function>

    INTEGER :: icontainer
    
    ! Which strategy should we use for searching for a memory block?
    SELECT CASE (rheap%ccontainerStrategy)
    CASE (EXSTOR_STRAT_LASTCONTAINER)
      ! Start at the end of the container list and proceed to the first container
      ! until we find one that can hold the data.
      DO icontainer = SIZE(rheap%p_RstorageContainers),1,-1
        IF ((rheap%p_RstorageContainers(icontainer)%imaxStorageMB .EQ. -1) .OR. &
            ((rheap%p_RstorageContainers(icontainer)%dcurrentStorage + &
              dsize)/1000000.0_DP .LT. &
              REAL(rheap%p_RstorageContainers(icontainer)%imaxStorageMB,DP))) THEN
          exstor_getContainer = icontainer
          RETURN
        END IF
      END DO
    
      ! Return the RamDrive container
      exstor_getContainer = 1
    END SELECT

  END FUNCTION

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_new1D (scall, sname, isize, ctype, ihandle, &
                           cinitNewBlock, icontainerId, rheap)

!<description>
  ! This routine reserves a 1D memory block of desired size and type.
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
  
  ! OPTIONAL: Id of the storage container (1,2,3,...) that should maintain 
  ! the memory block. If not specified, a container will automatically be 
  ! chosen (depending on the allocation strategy ccontainerStrategy of the
  ! container).
  ! Container-Id=1 always identifies the standard RamDrive container that
  ! saves data into the main memory. Container-Id=0 automatically determines
  ! a container (just as if the parameter is not specified).
  ! Note that if a valid container id <> 0 is specified, the memory block will
  ! be bound to that container and not be automatically moved to another one.
  INTEGER, INTENT(IN), OPTIONAL :: icontainerId

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  INTEGER, INTENT(OUT) :: ihandle

!</output>

!</subroutine>

    ! local variables
    INTEGER :: icontainer
    LOGICAL :: bbound
    REAL(DP) :: dmemsize
    CHARACTER(LEN=SYS_NAMELEN) :: snameBackup

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode

    IF (isize .EQ. 0) THEN
      CALL output_line ('isize=0', &
                        OU_CLASS_WARNING,OU_MODE_STD,'exstor_new1D')
      ihandle = ST_NOHANDLE
      RETURN
    END IF

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF
    
    ! Get the size of the memory block
    SELECT CASE (ctype)
    CASE (ST_SINGLE)
      dmemsize = REAL(isize,DP)*REAL(ST_SINGLE2BYTES)
    CASE (ST_DOUBLE)
      dmemsize = REAL(isize,DP)*REAL(ST_DOUBLE2BYTES)
    CASE (ST_INT)
      dmemsize = REAL(isize,DP)*REAL(ST_INT2BYTES)
    CASE (ST_LOGICAL)
      dmemsize = REAL(isize,DP)*REAL(ST_LOGICAL2BYTES)
    CASE (ST_CHAR)
      dmemsize = REAL(isize,DP)*REAL(ST_CHAR2BYTES)
    CASE DEFAULT
      CALL output_line ('Unknown mem type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_new1D')
      CALL sys_halt()
    END SELECT
    
    ! Determine the container id.
    icontainer = 0
    bbound = .FALSE.
    IF (PRESENT(icontainerId)) THEN
      icontainer = icontainerId
      
      ! The memory block is bound top that container if the
      ! container id is specified.
      bbound = icontainer .NE. 0
    END IF
    
    IF (icontainer .EQ. 0) THEN
      ! We have to find a container that is large enough to hold the data.
      ! Automatically determine a container.
      icontainer = exstor_getContainer (p_rheap,dmemsize)
    END IF
    
    IF ((icontainer .LT. 1) .OR. &
        (icontainer .GT. SIZE(p_rheap%p_RstorageContainers))) THEN
      CALL output_line ('Invalid container id!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_new1D')
      CALL sys_halt()
    END IF

    ! Back up the array name. This is important for a very crual situation:
    ! The storage_newhandle below may reallocate the p_Rdescriptors array.
    ! If sname points to one of the old arrays, the pointer gets invalid
    ! and the name cannot be accessed anymore. So make a backup of that 
    ! before creating a new handle!
    snameBackup = sname

    ! Get a new handle.
    ihandle = exstor_newhandle (p_rheap)

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Initialise the content

    p_rnode%cdataType = ctype
    p_rnode%idimension = 1
    p_rnode%sname = snameBackup
    p_rnode%icontainerId = icontainer
    p_rnode%bcontainerBound = bbound
    p_rnode%cinitNewBlock = cinitNewBlock
    p_rnode%dmemBytes = dmemsize
    p_rnode%Isize(1) = isize
    p_rnode%Isize(2) = 0

    p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
    IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem
      
    ! Notify the container about the new memory
    p_rheap%p_RstorageContainers(icontainer)%dcurrentStorage = &
        p_rheap%p_RstorageContainers(icontainer)%dcurrentStorage +  p_rnode%dmemBytes
        
    ! Some container-specific initialisation
    SELECT CASE (p_rheap%p_RstorageContainers(icontainer)%ctype)
    CASE (EXSTOR_CONT_DIRECTORY)
      ! Create a filename based on the filemane template of the container
      p_rnode%sfilename = &
          TRIM(p_rheap%p_RstorageContainers(icontainer)%sfilename)//'.'//&
          TRIM(sys_siL(ihandle,10))
    END SELECT

    ! Note: This routine will not immediately allocate memory!
    ! the actual memory allocation (or file creation) is done
    ! the first time, the memory is addressed!

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_new2D (scall, sname, Isize, ctype, ihandle, &
                           cinitNewBlock, icontainerId, rheap)

!<description>
  ! This routine reserves a 1D memory block of desired size and type.
!</description>

!<input>

  ! name of the calling routine
  CHARACTER(LEN=*), INTENT(IN) :: scall

  ! clear name of data field
  CHARACTER(LEN=*), INTENT(IN) :: sname

  ! requested storage size; DIMENSION(2)
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: Isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  INTEGER, INTENT(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT,ST_NEWBLOCK_ORDERED)
  INTEGER, INTENT(IN) :: cinitNewBlock
  
  ! OPTIONAL: Id of the storage container (1,2,3,...) that should maintain 
  ! the memory block. If not specified, a container will automatically be 
  ! chosen (depending on the allocation strategy ccontainerStrategy of the
  ! container).
  ! Container-Id=1 always identifies the standard RamDrive container that
  ! saves data into the main memory. Container-Id=0 automatically determines
  ! a container (just as if the parameter is not specified).
  ! Note that if a valid container id <> 0 is specified, the memory block will
  ! be bound to that container and not be automatically moved to another one.
  INTEGER, INTENT(IN), OPTIONAL :: icontainerId

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</inputoutput>

!<output>
  ! Handle of the memory block.
  INTEGER, INTENT(OUT) :: ihandle
!</output>

!</subroutine>

    ! local variables
    INTEGER :: icontainer
    LOGICAL :: bbound
    REAL(DP) :: dmemsize
    CHARACTER(LEN=SYS_NAMELEN) :: snameBackup

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode

    IF ((Isize(1) .EQ. 0) .OR. (Isize(2) .EQ. 0)) THEN
      CALL output_line ('Isize=0', &
                        OU_CLASS_WARNING,OU_MODE_STD,'exstor_new2D')
      ihandle = ST_NOHANDLE
      RETURN
    END IF

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF
    
    ! Get the size of the memory block
    SELECT CASE (ctype)
    CASE (ST_SINGLE)
      dmemsize = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_SINGLE2BYTES)
    CASE (ST_DOUBLE)
      dmemsize = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_DOUBLE2BYTES)
    CASE (ST_INT)
      dmemsize = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_INT2BYTES)
    CASE (ST_LOGICAL)
      dmemsize = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_LOGICAL2BYTES)
    CASE (ST_CHAR)
      dmemsize = REAL(Isize(1),DP)*REAL(Isize(2),DP)*REAL(ST_CHAR2BYTES)
    CASE DEFAULT
      CALL output_line ('Unknown mem type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_new2D')
      CALL sys_halt()
    END SELECT
    
    ! Determine the container id.
    icontainer = 0
    bbound = .FALSE.
    IF (PRESENT(icontainerId)) THEN
      icontainer = icontainerId
      
      ! The memory block is bound top that container if the
      ! container id is specified.
      bbound = icontainer .NE. 0
    END IF
    
    IF (icontainer .EQ. 0) THEN
      ! We have to find a container that is large enough to hold the data.
      ! Automatically determine a container.
      icontainer = exstor_getContainer (p_rheap,dmemsize)
    END IF
    
    IF ((icontainerId .LT. 1) .OR. &
        (icontainerId .GT. SIZE(p_rheap%p_RstorageContainers))) THEN
      CALL output_line ('Invalid container id!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_new2D')
      CALL sys_halt()
    END IF

    ! Back up the array name. This is important for a very crual situation:
    ! The storage_newhandle below may reallocate the p_Rdescriptors array.
    ! If sname points to one of the old arrays, the pointer gets invalid
    ! and the name cannot be accessed anymore. So make a backup of that 
    ! before creating a new handle!
    snameBackup = sname

    ! Get a new handle.
    ihandle = exstor_newhandle (p_rheap)

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Initialise the content

    p_rnode%cdataType = ctype
    p_rnode%idimension = 2
    p_rnode%sname = snameBackup
    p_rnode%icontainerId = icontainer
    p_rnode%bcontainerBound = bbound
    p_rnode%cinitNewBlock = cinitNewBlock
    p_rnode%dmemBytes = dmemsize
    p_rnode%Isize(1:2) = Isize(1:2)

    p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
    IF (p_rheap%dtotalMem .GT. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem
      
    ! Notify the container about the new memory
    p_rheap%p_RstorageContainers(icontainer)%dcurrentStorage = &
        p_rheap%p_RstorageContainers(icontainer)%dcurrentStorage +  p_rnode%dmemBytes

    ! Some container-specific initialisation
    SELECT CASE (p_rheap%p_RstorageContainers(icontainer)%ctype)
    CASE (EXSTOR_CONT_DIRECTORY)
      ! Create a filename based on the filemane template of the container
      p_rnode%sfilename = &
          TRIM(p_rheap%p_RstorageContainers(icontainer)%sfilename)//'.'//&
          TRIM(sys_siL(ihandle,10))
    END SELECT

    ! Note: This routine will not immediately allocate memory!
    ! the actual memory allocation (or file creation) is done
    ! the first time, the memory is addressed!

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_getsize1D (ihandle, isize, rheap)

!<description>
  ! Returns the length of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block 
  INTEGER, INTENT(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!<output>
  ! Length of the array identified by ihandle.
  INTEGER(I32), INTENT(OUT) :: isize
!</output>

!</subroutine>

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .LE. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getsize1D')
      CALL sys_halt()
    END IF

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Is the node associated at all?
    IF (p_rnode%cdataType .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid: '//sys_siL(ihandle,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getsize1D')
      CALL sys_halt()
    END IF

    ! What are we?
    IF (p_rnode%idimension .NE. 1) THEN
      CALL output_line ('Handle '//TRIM(sys_siL(ihandle,10))//' is not 1-dimensional!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getsize1D')
      CALL sys_halt()
    END IF
    
    isize = p_rnode%Isize(1)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_getsize2D (ihandle, isize, rheap)

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
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!<output>
  ! Length of each dimension of the array identified by ihandle.
  INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Isize
!</output>

!</subroutine>

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .LE. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      CALL sys_halt()
    END IF

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Is the node associated at all?
    IF (p_rnode%cdataType .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid: '//sys_siL(ihandle,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      CALL sys_halt()
    END IF

    ! What are we?
    IF (p_rnode%idimension .NE. 2) THEN
      CALL output_line ('Handle '//TRIM(sys_siL(ihandle,10))//' is not 2-dimensional!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      CALL sys_halt()
    END IF
    
    Isize(1:2) = p_rnode%Isize(1:2)
    
  END SUBROUTINE
  
!************************************************************************

!<subroutine>

  SUBROUTINE exstor_getdata_int (ihandle, Iarray, rheap)

!<description>
  ! This routine reads data from an external storage container and
  ! saves it to a local array.
!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! An array that is filled with the data identified by the handle.
  INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Iarray

!</output>

!</subroutine>

    ! local variables
    INTEGER :: iorder

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    IF (p_rheap%p_Rdescriptors(ihandle)%cdataType .NE. ST_INT) THEN
      CALL output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    IF (p_rnode%Isize(1) .NE. SIZE(Iarray)) THEN
      CALL output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    ! If the memory block is a pre-initialised block, we can directly
    ! fill it with data...    
    SELECT CASE (p_rnode%cinitNewBlock)
    CASE (ST_NEWBLOCK_ZERO)
      ! Fill the destination array with zero, that's all
      CALL lalg_clearVectorInt (Iarray)

    CASE (ST_NEWBLOCK_ORDERED)
      ! Fill the destination array by increasing numbers
      DO iorder=1,SIZE(Iarray)
        Iarray(iorder) = INT(iorder,I32)
      END DO

    CASE (ST_NEWBLOCK_NOINIT)
      ! Ok, there should be data behind. How to handle the data depends
      ! on the container type
      SELECT CASE(p_rcontainer%ctype)
      
      CASE (EXSTOR_CONT_RAMDRIVE)
        ! This is a RamDrive container. We use the memory management of the
        ! storage.f90 to maintain it.
        CALL getdata_ramdrive (p_rnode,Iarray)
                
      CASE (EXSTOR_CONT_DIRECTORY)
        ! This is a directory container maintaining the data as files on the 
        ! hard disc.
        CALL getdata_directory (p_rcontainer,p_rnode,Iarray)

      END SELECT
        
    END SELECT

  CONTAINS
  
    ! -------------------------------------------------------------------------
    SUBROUTINE getdata_ramdrive (rnode,dataarray)
    
    ! Retrieves the data from a RamDrive container
    
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(IN) :: rnode
    
    ! The destination array for the data
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: dataarray
    
      ! local variables
      INTEGER(I32), DIMENSION(:), POINTER :: p_data

      ! At first: Do we have data at all?
      IF (rnode%istorageHandle .EQ. ST_NOHANDLE) THEN
        CALL output_line ('Trying to read noninitialised memory! Handle: '//&
            sys_siL(ihandle,10), OU_CLASS_ERROR,OU_MODE_STD,'getdata_ramdrive')
        CALL sys_halt()
      END IF
    
      ! Get the memory block and copy it.
      CALL storage_getbase_int (rnode%istoragehandle,p_data)
      CALL lalg_copyVectorInt (p_data,dataarray)
      
    END SUBROUTINE

    ! -------------------------------------------------------------------------
    SUBROUTINE getdata_directory (rcontainer,rnode,dataarray)
    
    ! Retrieves the data from a directory container
    
    ! The storage container assigned to the storage block
    TYPE(t_exstorageContainer), INTENT(IN) :: rcontainer
        
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(IN) :: rnode
    
    ! The destination array for the data
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: dataarray
    
      ! local variables
      INTEGER :: cf
      
      ! Open the file to read data
      CALL io_openFileForReading(TRIM(rcontainer%spath)//rnode%sfilename, &
          cf, bformatted=rcontainer%bformatted)
      
      ! Read the data from the file
      IF (rcontainer%bformatted) THEN
        READ(cf,*) dataarray(:)
      ELSE
        READ(cf) dataarray(:)
      END IF
      
      ! Close the file, finish
      CLOSE (cf)
      
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_setdata_int (ihandle, Iarray, rheap)

!<description>
  ! This routine writes data from a local array top an external storage 
  ! container.
!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

  ! An array with data that should be stored on the external storage
  ! container.
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: Iarray
!</input>

!</subroutine>

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    IF (p_rheap%p_Rdescriptors(ihandle)%cdataType .NE. ST_INT) THEN
      CALL output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    IF (p_rnode%Isize(1) .NE. SIZE(Iarray)) THEN
      CALL output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    ! The data block will now contain data.
    p_rnode%cinitNewBlock = ST_NEWBLOCK_NOINIT

    ! How to handle the data depends on the container type.
    SELECT CASE(p_rcontainer%ctype)
    
    CASE (EXSTOR_CONT_RAMDRIVE)
      ! This is a RamDrive container. We use the memory management of the
      ! storage.f90 to maintain it.
      CALL setdata_ramdrive (p_rnode,Iarray)
              
    CASE (EXSTOR_CONT_DIRECTORY)
      ! This is a directory container maintaining the data as files on the 
      ! hard disc.
      CALL setdata_directory (p_rcontainer,p_rnode,Iarray)

    END SELECT
        
  CONTAINS
  
    ! -------------------------------------------------------------------------
    SUBROUTINE setdata_ramdrive (rnode,dataarray)
    
    ! Writes data to a RamDrive container
    
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(INOUT) :: rnode
    
    ! The source array with the data
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: dataarray
    
      ! local variables
      INTEGER(I32), DIMENSION(:), POINTER :: p_data

      ! Copy data with the storage-copy command. If necessary, new memory
      ! is allocated.
      IF (rnode%istorageHandle .EQ. ST_NOHANDLE) THEN
        ! We have to allocate memory.
        CALL storage_new ('setdata_ramdrive','exdata',rnode%Isize(1),&
            rnode%cdataType,rnode%istorageHandle,ST_NEWBLOCK_NOINIT)
      END IF
    
      ! Get the memory block and copy it.
      CALL storage_getbase_int (rnode%istoragehandle,p_data)
      CALL lalg_copyVectorInt (dataarray,p_data)
      
    END SUBROUTINE

    ! -------------------------------------------------------------------------
    SUBROUTINE setdata_directory (rcontainer,rnode,dataarray)
    
    ! Writes data to a directory container
    
    ! The storage container assigned to the storage block
    TYPE(t_exstorageContainer), INTENT(IN) :: rcontainer
        
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(INOUT) :: rnode
    
    ! The source array with the data
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: dataarray
    
      ! local variables
      INTEGER :: cf
      
      ! Open the file to write data
      CALL io_openFileForWriting(TRIM(rcontainer%spath)//rnode%sfilename, &
          cf, SYS_REPLACE, bformatted=rcontainer%bformatted)
      
      ! Write the data to the file
      IF (rcontainer%bformatted) THEN
        WRITE(cf,*) dataarray(:)
      ELSE
        WRITE(cf) dataarray(:) 
      END IF
      
      ! Close the file, finish
      CLOSE (cf)
      
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_getdata_single (ihandle, Farray, rheap)

!<description>
  ! This routine reads data from an external storage container and
  ! saves it to a local array.
!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! An array that is filled with the data identified by the handle.
  REAL(SP), DIMENSION(:), INTENT(OUT) :: Farray

!</output>

!</subroutine>

    ! local variables
    INTEGER :: iorder

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    IF (p_rheap%p_Rdescriptors(ihandle)%cdataType .NE. ST_SINGLE) THEN
      CALL output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    IF (p_rnode%Isize(1) .NE. SIZE(Farray)) THEN
      CALL output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    ! If the memory block is a pre-initialised block, we can directly
    ! fill it with data...    
    SELECT CASE (p_rnode%cinitNewBlock)
    CASE (ST_NEWBLOCK_ZERO)
      ! Fill the destination array with zero, that's all
      CALL lalg_clearVectorSngl (Farray)

    CASE (ST_NEWBLOCK_ORDERED)
      ! Fill the destination array by increasing numbers
      DO iorder=1,SIZE(Farray)
        Farray(iorder) = REAL(iorder,SP)
      END DO

    CASE (ST_NEWBLOCK_NOINIT)
      ! Ok, there should be data behind. How to handle the data depends
      ! on the container type
      SELECT CASE(p_rcontainer%ctype)
      
      CASE (EXSTOR_CONT_RAMDRIVE)
        ! This is a RamDrive container. We use the memory management of the
        ! storage.f90 to maintain it.
        CALL getdata_ramdrive (p_rnode,Farray)
                
      CASE (EXSTOR_CONT_DIRECTORY)
        ! This is a directory container maintaining the data as files on the 
        ! hard disc.
        CALL getdata_directory (p_rcontainer,p_rnode,Farray)

      END SELECT
        
    END SELECT

  CONTAINS
  
    ! -------------------------------------------------------------------------
    SUBROUTINE getdata_ramdrive (rnode,dataarray)
    
    ! Retrieves the data from a RamDrive container
    
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(IN) :: rnode
    
    ! The destination array for the data
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dataarray
    
      ! local variables
      REAL(SP), DIMENSION(:), POINTER :: p_data

      ! At first: Do we have data at all?
      IF (rnode%istorageHandle .EQ. ST_NOHANDLE) THEN
        CALL output_line ('Trying to read noninitialised memory! Handle: '//&
            sys_siL(ihandle,10), OU_CLASS_ERROR,OU_MODE_STD,'getdata_ramdrive')
        CALL sys_halt()
      END IF
    
      ! Get the memory block and copy it.
      CALL storage_getbase_single (rnode%istoragehandle,p_data)
      CALL lalg_copyVectorSngl (p_data,dataarray)
      
    END SUBROUTINE

    ! -------------------------------------------------------------------------
    SUBROUTINE getdata_directory (rcontainer,rnode,dataarray)
    
    ! Retrieves the data from a directory container
    
    ! The storage container assigned to the storage block
    TYPE(t_exstorageContainer), INTENT(IN) :: rcontainer
        
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(IN) :: rnode
    
    ! The destination array for the data
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dataarray
    
      ! local variables
      INTEGER :: cf
      
      ! Open the file to read data
      CALL io_openFileForReading(TRIM(rcontainer%spath)//rnode%sfilename, &
          cf, bformatted=rcontainer%bformatted)
      
      ! Read the data from the file
      IF (rcontainer%bformatted) THEN
        READ(cf,*) dataarray(:)
      ELSE
        READ(cf) dataarray(:)
      END IF
      
      ! Close the file, finish
      CLOSE (cf)
      
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_setdata_single (ihandle, Farray, rheap)

!<description>
  ! This routine writes data from a local array top an external storage 
  ! container.
!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

  ! An array with data that should be stored on the external storage
  ! container.
  REAL(SP), DIMENSION(:), INTENT(IN) :: Farray
!</input>

!</subroutine>

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    IF (p_rheap%p_Rdescriptors(ihandle)%cdataType .NE. ST_SINGLE) THEN
      CALL output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    IF (p_rnode%Isize(1) .NE. SIZE(Farray)) THEN
      CALL output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    ! The data block will now contain data.
    p_rnode%cinitNewBlock = ST_NEWBLOCK_NOINIT

    ! How to handle the data depends on the container type.
    SELECT CASE(p_rcontainer%ctype)
    
    CASE (EXSTOR_CONT_RAMDRIVE)
      ! This is a RamDrive container. We use the memory management of the
      ! storage.f90 to maintain it.
      CALL setdata_ramdrive (p_rnode,Farray)
              
    CASE (EXSTOR_CONT_DIRECTORY)
      ! This is a directory container maintaining the data as files on the 
      ! hard disc.
      CALL setdata_directory (p_rcontainer,p_rnode,Farray)

    END SELECT
        
  CONTAINS
  
    ! -------------------------------------------------------------------------
    SUBROUTINE setdata_ramdrive (rnode,dataarray)
    
    ! Writes data to a RamDrive container
    
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(INOUT) :: rnode
    
    ! The source array with the data
    REAL(SP), DIMENSION(:), INTENT(IN) :: dataarray
    
      ! local variables
      REAL(SP), DIMENSION(:), POINTER :: p_data

      ! Copy data with the storage-copy command. If necessary, new memory
      ! is allocated.
      IF (rnode%istorageHandle .EQ. ST_NOHANDLE) THEN
        ! We have to allocate memory.
        CALL storage_new ('setdata_ramdrive','exdata',rnode%Isize(1),&
            rnode%cdataType,rnode%istorageHandle,ST_NEWBLOCK_NOINIT)
      END IF
    
      ! Get the memory block and copy it.
      CALL storage_getbase_single (rnode%istoragehandle,p_data)
      CALL lalg_copyVectorSngl (dataarray,p_data)
      
    END SUBROUTINE

    ! -------------------------------------------------------------------------
    SUBROUTINE setdata_directory (rcontainer,rnode,dataarray)
    
    ! Writes data to a directory container
    
    ! The storage container assigned to the storage block
    TYPE(t_exstorageContainer), INTENT(IN) :: rcontainer
        
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(INOUT) :: rnode
    
    ! The source array with the data
    REAL(SP), DIMENSION(:), INTENT(IN) :: dataarray
    
      ! local variables
      INTEGER :: cf
      
      ! Open the file to write data
      CALL io_openFileForWriting(TRIM(rcontainer%spath)//rnode%sfilename, &
          cf, SYS_REPLACE, bformatted=rcontainer%bformatted)
      
      ! Write the data to the file
      IF (rcontainer%bformatted) THEN
        WRITE(cf,*) dataarray(:)
      ELSE
        WRITE(cf) dataarray(:) 
      END IF
      
      ! Close the file, finish
      CLOSE (cf)
      
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_getdata_double (ihandle, Darray, rheap)

!<description>
  ! This routine reads data from an external storage container and
  ! saves it to a local array.
!</description>

!<input>

  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

!</input>

!<output>

  ! An array that is filled with the data identified by the handle.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Darray

!</output>

!</subroutine>

    ! local variables
    INTEGER :: iorder

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    IF (p_rheap%p_Rdescriptors(ihandle)%cdataType .NE. ST_DOUBLE) THEN
      CALL output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    IF (p_rnode%Isize(1) .NE. SIZE(Darray)) THEN
      CALL output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      CALL sys_halt()
    END IF

    ! If the memory block is a pre-initialised block, we can directly
    ! fill it with data...    
    SELECT CASE (p_rnode%cinitNewBlock)
    CASE (ST_NEWBLOCK_ZERO)
      ! Fill the destination array with zero, that's all
      CALL lalg_clearVectorDble (Darray)

    CASE (ST_NEWBLOCK_ORDERED)
      ! Fill the destination array by increasing numbers
      DO iorder=1,SIZE(Darray)
        Darray(iorder) = REAL(iorder,DP)
      END DO

    CASE (ST_NEWBLOCK_NOINIT)
      ! Ok, there should be data behind. How to handle the data depends
      ! on the container type
      SELECT CASE(p_rcontainer%ctype)
      
      CASE (EXSTOR_CONT_RAMDRIVE)
        ! This is a RamDrive container. We use the memory management of the
        ! storage.f90 to maintain it.
        CALL getdata_ramdrive (p_rnode,Darray)
                
      CASE (EXSTOR_CONT_DIRECTORY)
        ! This is a directory container maintaining the data as files on the 
        ! hard disc.
        CALL getdata_directory (p_rcontainer,p_rnode,Darray)

      END SELECT
        
    END SELECT

  CONTAINS
  
    ! -------------------------------------------------------------------------
    SUBROUTINE getdata_ramdrive (rnode,dataarray)
    
    ! Retrieves the data from a RamDrive container
    
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(IN) :: rnode
    
    ! The destination array for the data
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dataarray
    
      ! local variables
      REAL(DP), DIMENSION(:), POINTER :: p_data

      ! At first: Do we have data at all?
      IF (rnode%istorageHandle .EQ. ST_NOHANDLE) THEN
        CALL output_line ('Trying to read noninitialised memory! Handle: '//&
            sys_siL(ihandle,10), OU_CLASS_ERROR,OU_MODE_STD,'getdata_ramdrive')
        CALL sys_halt()
      END IF
    
      ! Get the memory block and copy it.
      CALL storage_getbase_double (rnode%istoragehandle,p_data)
      CALL lalg_copyVectorDble (p_data,dataarray)
      
    END SUBROUTINE

    ! -------------------------------------------------------------------------
    SUBROUTINE getdata_directory (rcontainer,rnode,dataarray)
    
    ! Retrieves the data from a directory container
    
    ! The storage container assigned to the storage block
    TYPE(t_exstorageContainer), INTENT(IN) :: rcontainer

    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(IN) :: rnode
    
    ! The destination array for the data
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dataarray
    
      ! local variables
      INTEGER :: cf
      
      ! Open the file to read data
      CALL io_openFileForReading(TRIM(rcontainer%spath)//rnode%sfilename, &
          cf, bformatted=rcontainer%bformatted)
      
      ! Read the data from the file
      IF (rcontainer%bformatted) THEN
        READ(cf,*) dataarray(:)
      ELSE
        READ(cf) dataarray(:)
      END IF
      
      ! Close the file, finish
      CLOSE (cf)
      
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_setdata_double (ihandle, Darray, rheap)

!<description>
  ! This routine writes data from a local array top an external storage 
  ! container.
!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

  ! An array with data that should be stored on the external storage
  ! container.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Darray
!</input>

!</subroutine>

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    IF (p_rheap%p_Rdescriptors(ihandle)%cdataType .NE. ST_DOUBLE) THEN
      CALL output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    IF (p_rnode%Isize(1) .NE. SIZE(Darray)) THEN
      CALL output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      CALL sys_halt()
    END IF

    ! The data block will now contain data.
    p_rnode%cinitNewBlock = ST_NEWBLOCK_NOINIT

    ! How to handle the data depends on the container type.
    SELECT CASE(p_rcontainer%ctype)
    
    CASE (EXSTOR_CONT_RAMDRIVE)
      ! This is a RamDrive container. We use the memory management of the
      ! storage.f90 to maintain it.
      CALL setdata_ramdrive (p_rnode,Darray)
              
    CASE (EXSTOR_CONT_DIRECTORY)
      ! This is a directory container maintaining the data as files on the 
      ! hard disc.
      CALL setdata_directory (p_rcontainer,p_rnode,Darray)

    END SELECT
        
  CONTAINS
  
    ! -------------------------------------------------------------------------
    SUBROUTINE setdata_ramdrive (rnode,dataarray)
    
    ! Writes data to a RamDrive container
    
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(INOUT) :: rnode
    
    ! The source array with the data
    REAL(DP), DIMENSION(:), INTENT(IN) :: dataarray
    
      ! local variables
      REAL(DP), DIMENSION(:), POINTER :: p_data

      ! Copy data with the storage-copy command. If necessary, new memory
      ! is allocated.
      IF (rnode%istorageHandle .EQ. ST_NOHANDLE) THEN
        ! We have to allocate memory.
        CALL storage_new ('setdata_ramdrive','exdata',rnode%Isize(1),&
            rnode%cdataType,rnode%istorageHandle,ST_NEWBLOCK_NOINIT)
      END IF
    
      ! Get the memory block and copy it.
      CALL storage_getbase_double (rnode%istoragehandle,p_data)
      CALL lalg_copyVectorDble (dataarray,p_data)
      
    END SUBROUTINE

    ! -------------------------------------------------------------------------
    SUBROUTINE setdata_directory (rcontainer,rnode,dataarray)
    
    ! Writes data to a directory container
    
    ! The storage container assigned to the storage block
    TYPE(t_exstorageContainer), INTENT(IN) :: rcontainer
        
    ! The storage block containing the data
    TYPE(t_exstorageNode), INTENT(INOUT) :: rnode
    
    ! The source array with the data
    REAL(DP), DIMENSION(:), INTENT(IN) :: dataarray
    
      ! local variables
      INTEGER :: cf
      
      ! Open the file to write data
      CALL io_openFileForWriting(TRIM(rcontainer%spath)//rnode%sfilename, &
          cf, SYS_REPLACE, bformatted=rcontainer%bformatted)
      
      ! Write the data to the file
      IF (rcontainer%bformatted) THEN
        WRITE(cf,*) dataarray(:)
      ELSE
        WRITE(cf) dataarray(:) 
      END IF
      
      ! Close the file, finish
      CLOSE (cf)
      
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_getdata_storage (ihandle, istoragehandle, rheap, rstorageheap)

!<description>
  ! This routine reads data from an external storage container and
  ! saves it to a local array in the storage management of the main
  ! memory.
!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle
  
  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</input>

!<inputoutput>
  ! Handle to a memory block in the main memory, maintained by the 
  ! storage.f90. Data read from the data container is directly saved
  ! to this memory block.
  ! If this is ST_NOHANDLE, a new handle is automatically created
  ! and the data is saved to it.
  INTEGER, INTENT(INOUT) :: istoragehandle

  ! OPTIONAL: local heap structure of the storage management for the
  ! main memory. If not given, the global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rstorageheap
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ctype
    INTEGER(I32) :: isize

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
      CALL sys_halt()
    END IF

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    IF (istoragehandle .EQ. ST_NOHANDLE) THEN
    
      ! Allocate a new memory block with the correct data type that is just
      ! as large that it can hold our data.
      
      SELECT CASE(p_rnode%idimension)
      CASE (1)
        CALL storage_new ('exstor_getdata_storage','datacopy',p_rnode%Isize(1),&
            p_rnode%cdataType,istoragehandle,ST_NEWBLOCK_NOINIT,rstorageheap)
      CASE (2)
        CALL storage_new ('exstor_getdata_storage','datacopy',p_rnode%Isize,&
            p_rnode%cdataType,istoragehandle,ST_NEWBLOCK_NOINIT,rstorageheap)
      END SELECT
      
    ELSE
      CALL storage_getdatatype(istoragehandle,ctype)
      IF (p_rnode%cdataType .NE. ctype) THEN
        CALL output_line ('Wrong data format!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
        CALL sys_halt()
      END IF
      
      CALL storage_getsize (istoragehandle,isize)
      IF (p_rnode%Isize(1) .NE. isize) THEN
        CALL output_line ('Data array has the wrong size!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
        CALL sys_halt()
      END IF

    END IF
    
    ! If the memory block is a pre-initialised block, we can directly
    ! fill it with data...    
    IF (p_rnode%cinitNewBlock .NE. ST_NEWBLOCK_NOINIT) THEN
      CALL storage_initialiseBlock (istoragehandle, p_rnode%cinitNewBlock)
      RETURN
    END IF
    
    ! Ok, we have to copy data from the external storage to the memory block.
    ! This now depends on the data type...
    SELECT CASE (p_rnode%cdataType)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double (istoragehandle,p_Ddata,rstorageheap)
      CALL exstor_getdata_double (ihandle,p_Ddata)
    CASE (ST_SINGLE)
      CALL storage_getbase_single (istoragehandle,p_Fdata,rstorageheap)
      CALL exstor_getdata_single (ihandle,p_Fdata)
    CASE (ST_INT)
      CALL storage_getbase_int (istoragehandle,p_Idata,rstorageheap)
      CALL exstor_getdata_int (ihandle,p_Idata)
    CASE DEFAULT
      CALL output_line ('Unsupported data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
    END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_setdata_storage (ihandle, istoragehandle, rheap, rstorageheap)

!<description>
  ! This routine writes data from a local array in the storage management
  ! of the main memory to an external storage container.
!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle
  
  ! Handle to a memory block in the main memory, maintained by the 
  ! storage.f90. The data in this memory block is directly saved
  ! to the the data container.
  INTEGER, INTENT(IN) :: istoragehandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap

  ! OPTIONAL: local heap structure of the storage management for the
  ! main memory. If not given, the global heap is used.
  TYPE(t_storageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rstorageheap
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ctype
    INTEGER(I32) :: isize

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    TYPE(t_exstorageContainer), POINTER :: p_rcontainer
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
      CALL sys_halt()
    END IF

    IF (istoragehandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('istoragehandle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
      CALL sys_halt()
    END IF

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    CALL storage_getdatatype(istoragehandle,ctype)
    IF (p_rnode%cdataType .NE. ctype) THEN
      CALL output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
      CALL sys_halt()
    END IF
    
    CALL storage_getsize (istoragehandle,isize)
    IF (p_rnode%Isize(1) .NE. isize) THEN
      CALL output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
      CALL sys_halt()
    END IF
    
    ! Ok, we have to copy data from the external storage to the memory block.
    ! This now depends on the data type...
    SELECT CASE (p_rnode%cdataType)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double (istoragehandle,p_Ddata)
      CALL exstor_setdata_double (ihandle,p_Ddata)
    CASE (ST_SINGLE)
      CALL storage_getbase_single (istoragehandle,p_Fdata)
      CALL exstor_setdata_single (ihandle,p_Fdata)
    CASE (ST_INT)
      CALL storage_getbase_int (istoragehandle,p_Idata)
      CALL exstor_setdata_int (ihandle,p_Idata)
    CASE DEFAULT
      CALL output_line ('Unsupported data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
    END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_clear (ihandle, rheap)

!<description>
  ! This routine clears an array identified by ihandle; all entries are
  ! overwritten by 0.
!</description>

!<input>
  ! The handle
  INTEGER, INTENT(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
!</inputoutput>

!</subroutine>

    ! Pointer to the heap
    TYPE(t_exstorageBlock), POINTER :: p_rheap

    ! Get the heap to use - local or global one.

    IF(PRESENT(rheap)) THEN
      p_rheap => rheap
    ELSE
      p_rheap => rbaseexternal
    END IF

    IF (ihandle .EQ. ST_NOHANDLE) THEN
      CALL output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_clear')
      CALL sys_halt()
    END IF

    ! Set the initialised-flag to ST_NEWBLOCK_ZERO. The data is
    ! trated as zero when it's accessed the next time.
    p_rheap%p_Rdescriptors(ihandle)%cinitNewBlock = ST_NEWBLOCK_ZERO

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE exstor_copy(h_source, h_dest, icontainerId, rheap)

!<description>
  ! This routine copies the information of one array to another.
  ! The structure of the arrays behind h_source and h_dest must be the
  ! same!
!</description>

!<input>
  ! Handle of the source array to copy
  INTEGER, INTENT(IN) :: h_source

  ! OPTIONAL: If h_dest=ST_NOHANDLE, this parameter allows to specify
  ! a container id where new data is stored to.
  ! If h_dest<>ST_NOHANDLE, this parameter is ignored.
  INTEGER, INTENT(IN), OPTIONAL :: icontainerId

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  TYPE(t_exstorageBlock), INTENT(INOUT), TARGET, OPTIONAL :: rheap
  
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it. Such a new
  ! block is allocated in the container icontainerId
  INTEGER, INTENT(INOUT) :: h_dest
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_exstorageBlock), POINTER :: p_rheap
    TYPE(t_exstorageNode), POINTER :: p_rnode
    INTEGER :: ihandle
    
    ! Read in the data, create a temp array in the main memory.
    ihandle = ST_NOHANDLE
    CALL exstor_getdata_storage (h_source, ihandle, rheap)
    
    ! If necessary, allocate a new block.
    IF (h_dest .EQ. ST_NOHANDLE) THEN

      IF(PRESENT(rheap)) THEN
        p_rheap => rheap
      ELSE
        p_rheap => rbaseexternal
      END IF

      IF (h_source .EQ. ST_NOHANDLE) THEN
        CALL output_line ('Handle invalid!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
        CALL sys_halt()
      END IF
    
      ! Get the container and the storage node
      p_rnode => p_rheap%p_Rdescriptors(ihandle)
      
      ! Allocate a new memory block with the correct data type that is just
      ! as large that it can hold our data.
      
      SELECT CASE(p_rnode%idimension)
      CASE (1)
        CALL exstor_new ('exstor_copy','datacopy',p_rnode%Isize(1),&
            p_rnode%cdataType,h_dest,ST_NEWBLOCK_NOINIT,icontainerId,rheap)
      CASE (2)
        CALL exstor_new ('exstor_copy','datacopy',p_rnode%Isize,&
            p_rnode%cdataType,h_dest,ST_NEWBLOCK_NOINIT,icontainerId,rheap)
      END SELECT
      
    END IF
    
    ! Write the data to a new block
    CALL exstor_setdata_storage (h_dest, ihandle, rheap)
    
    ! Release the temp memory
    CALL storage_free (ihandle)

  END SUBROUTINE

END MODULE
