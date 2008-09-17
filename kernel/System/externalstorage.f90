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

module externalstorage

  use fsystem
  use genoutput
  use storage
  use io

  implicit none

!<constants>

!<constantblock description="Type flags identifying external storage container.">

  ! Undefined container
  integer, parameter :: EXSTOR_CONT_UNDEFINED   = -1

  ! In-memory/RamDrive container. The 'external storage' is actually 'in memory' and
  ! not on any hard disc or similar external storage.
  integer, parameter :: EXSTOR_CONT_RAMDRIVE    = 0
  
  ! The external storage container is a subdirectory on the hard disc.
  integer, parameter :: EXSTOR_CONT_DIRECTORY   = 1

!</constantblock>

!<constantblock description="Strategy constants for choosing a container.">

  ! Always choose the last existing container. If it's full, choose the
  ! last but one, etc.
  integer, parameter :: EXSTOR_STRAT_LASTCONTAINER = 0

!</constantblock>

!</constants>

!<types>

  !<typeblock>

  ! This type block specifies a container for external storage.
  ! This may be for example a directory on a hard disc.

  type t_exstorageContainer

    private
    
    ! Type of the external storage container. One of the EXSTOR_CONT_xxxx
    ! flags.
    integer :: ctype = EXSTOR_CONT_UNDEFINED
    
    ! Maximum storage (in Megabytes) that this storage container can handle.
    ! =-1: Infinite memory available.
    integer :: imaxStorageMB = -1
    
    ! Number of bytes currently allocated by this storage container.
    ! This is an integer number encoded as a double precision number
    ! in order to capture even the size of large memory drives.
    real(DP) :: dcurrentStorage = 0.0_DP
    
    ! For ctype=EXSTOR_CONT_DIRECTORY: Name of the directory containing
    ! the files with data. The string contains a "/" at the end.
    character(LEN=SYS_STRLEN) :: spath = './'
    
    ! For ctype=EXSTOR_CONT_DIRECTORY: Basic filename of the data files.
    character(LEN=SYS_NAMELEN) :: sfilename = 'feat2exstorage'
    
    ! For ctype=EXSTOR_CONT_DIRECTORY: Whether to save the data formatted 
    ! or unformatted.
    logical :: bformatted = .false.

  end type

  !</typeblock>

  !<typeblock>

  ! Type block for describing a handle. This collects the number of the
  ! handle, the storage amount associated to it, the pointer to the memory
  ! location etc.

  type t_exstorageNode
  
    private

    ! Type of data associated to the handle (ST_NOHANDLE, ST_SINGLE,
    ! ST_DOUBLE, ST_INT, ST_LOGICAL, ST_CHAR)
    integer :: cdataType = ST_NOHANDLE

    ! Dimension associated to the handle (0=not assigned, 1=1D, 2=2D array)
    integer :: idimension = 0

    ! The name of the array that is associated to that handle
    character(LEN=SYS_NAMELEN) :: sname = ''

    ! Amount of memory (in bytes) associated to this block.
    ! We store that as a double to allow storing numbers > 2GB !
    real(DP) :: dmemBytes = 0.0_DP
    
    ! Size of the data array. For 1D arrays, the array size can be found
    ! in Isize(1).
    integer(I32), dimension(2) :: Isize = (/0,0/)

    ! Id of the external storage container that contains the data
    integer :: icontainerId = 0
    
    ! Whether the storage block is bound to a container or not.
    ! TRUE=The memory block is bound to container icontainerId.
    ! FALSE=The memory management system may automatically move the memory 
    !       block from one container to another if necessary.
    logical :: bcontainerBound = .false.

    ! Flag that identifies whether this storage block is initialised
    ! somehow.
    ! 
    ! Flag that specifies whether the data in this vector is 
    ! treated as being completely zero.
    ! ST_NEWBLOCK_NOINIT:  The storage block contains data.
    ! ST_NEWBLOCK_ZERO:    The storage block is to be treated as zero.
    ! ST_NEWBLOCK_ORDERED: The storage block is to be treated as initialised
    !                      by a sequence of numbers (1,2,3,...)
    integer :: cinitNewBlock = ST_NEWBLOCK_NOINIT
    
    ! This is a handle to a memory block. The handle is valid if the
    ! storage block is realised in the main memory as part of a RamDrive
    ! container. (A RamDrive container uses the standard memory management 
    ! of the storage.f90 and saves memory blocks in the main memory.)
    integer :: istorageHandle = ST_NOHANDLE
    
    ! if the storage block is realised as a file in a directory, this
    ! variable contains the filename of the file.
    character(LEN=SYS_NAMELEN) :: sfilename = ''
    
  end type

  !</typeblock>

  !<typeblock>

  ! This block represents a heap that maintains singole, double precision
  ! and integer data. It contains a list of t_storageNode elements for all
  ! the handles.
  ! There's one global object of this type for the global storage management,
  ! but if necessary, an algorithm can create such a block locally, too,
  ! to prevent conflicts with the global memory.

  type t_exstorageBlock

    private

    ! An array of t_exstorageContainer objects that identifies all
    ! possible storage containers that are handled by this storage
    ! block.
    type(t_exstorageContainer), dimension(:), pointer :: p_RstorageContainers => null()

    ! Strategy how to choose a container in p_RstorageContainers where
    ! to save data. One of the EXSTOR_STRAT_xxxx costants. By default,
    ! this is EXSTOR_STRAT_LASTCONTAINER, so always the last container
    ! providing enough free memory is chosen to save new data.
    integer :: ccontainerStrategy = EXSTOR_STRAT_LASTCONTAINER

    ! An array of t_exstorageNode objects corresponding to the handles.
    ! Each entry identifies a memory block in an external storage container.
    ! Can be dynamically extended if there are not enough handles available.
    type(t_exstorageNode), dimension(:), pointer :: p_Rdescriptors => null()

    ! A list of all 'free' handles. This is a 'ring' queue. If all
    ! handles are in use, p_Rdescriptors and p_IfreeHandles are dynamically
    ! extended.
    integer, dimension(:), pointer :: p_IfreeHandles => null()

    ! Index in p_IfreeHandles to the next free handle
    integer :: p_inextFreeHandle = 0

    ! Index in p_IfreeHandles to the last free handle
    integer :: p_ilastFreeHandle = 0

    ! Number of handles in use
    integer :: ihandlesInUse = 0

    ! Total number of handles maintained by this block; = size(p_Rdescriptors).
    integer :: nhandlesTotal = 0

    ! Number of handles to add if there are not enough free handles.
    integer :: ihandlesDelta = 0

    ! Total amount of memory (in bytes) that is in use. We maintain it
    ! as a double as this allows to save values > 2GB!
    real(DP) :: dtotalMem = 0.0_DP

    ! Maximum number of handles that were in use ofer the whole lifetime
    ! of this structure.
    integer :: nhandlesInUseMax = 0

    ! Maximum amount of memory that was in use over the whole lifetime
    ! of this structure.
    real(DP) :: dtotalMemMax = 0.0_DP

  end type

  !</typeblock>

!</types>

!<globals>

  ! Global memory management structure
  type(t_exstorageBlock), save, target :: rbaseExternal

!</globals>

  interface exstor_new
    module procedure exstor_new1D
    module procedure exstor_new2D
  end interface
  
  interface exstor_getsize
    module procedure exstor_getsize1D
    module procedure exstor_getsize2D
  end interface
  
  private :: exstor_newhandle
  private :: exstor_releasehandle
  private :: exstor_getContainer

contains

!************************************************************************

!<subroutine>

  subroutine exstor_init(ihandleCount, ihandlesDelta, rheap)

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
  integer, intent(IN) :: ihandleCount

  ! OPTIONAL: Number of handles to increase the memory block by, if there are
  ! not enough handles available. Standard setting is 1/2*ihandleCount.
  integer, intent(IN), optional :: ihandlesDelta

!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is initialised.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!</subroutine>

    ! local variables

    ! the real 'handle-delta'
    integer :: ihandles, ihDelta

    ! Pointer to the heap to initialise
    type(t_exstorageBlock), pointer :: p_rheap

    integer :: i

    ! Initialise ihDelta and p_rheap and work with these - as the other
    ! parameters are optional.
    ! We work at least with 1 handles and ihDelta = 1.

    ihandles = max(1,ihandlecount)

    ihDelta = 1
    if(present(ihandlesDelta)) ihDelta = ihandlesDelta
    ihDelta = max(1,ihDelta)

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    ! Initialise the memory management block

    p_rheap%nhandlesTotal = ihandles
    p_rheap%ihandlesDelta = ihDelta
    p_rheap%p_inextFreeHandle = 1
    p_rheap%p_ilastFreeHandle = ihandles
    p_rheap%ihandlesInUse = 0
    p_rheap%nhandlesInUseMax = 0
    allocate(p_rheap%p_Rdescriptors(ihandles))
    allocate(p_rheap%p_IfreeHandles(ihandles))

    ! All handles free
    do i=1,ihandles
      p_rheap%p_IfreeHandles(i) = i
    end do
    
    ! Attach a RamDrive container with arbitrary memory size.
    call exstor_attachRamdrive(-1,p_rheap)

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_attachRamdrive(imaxMem,rheap)

!<description>
  ! This routine attaches a Ramdrive container to the external storage
  ! management.
!</description>

!<input>
  ! OPTIONAL: Maximum size of the RamDrive (in Megabytes).
  ! -1 or not present = arbitrary.
  integer, intent(IN), optional :: imaxMem
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is initialised.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_exstorageContainer), dimension(:), pointer :: p_RstorageContainers
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Pointer to the heap to initialise
    type(t_exstorageBlock), pointer :: p_rheap

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    ! Create a new container.
    if (.not. associated(p_rheap%p_RstorageContainers)) then
      allocate(p_RstorageContainers(1))
    else
      allocate(p_RstorageContainers(size(p_rheap%p_RstorageContainers)+1))
      p_RstorageContainers(1:size(p_rheap%p_RstorageContainers)) = &
          p_rheap%p_RstorageContainers(:)
      
      deallocate(p_rheap%p_RstorageContainers)
    end if
    p_rheap%p_RstorageContainers => p_RstorageContainers
    
    p_rcontainer => p_rheap%p_RstorageContainers(size(p_rheap%p_RstorageContainers))
    
    ! Initialise the container as a RamDrive container
    p_rcontainer%ctype = EXSTOR_CONT_RAMDRIVE
    if(present(imaxMem)) p_rcontainer%imaxStorageMB = imaxMem
    
  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_attachDirectory(spath,imaxMem,sfilename,bformatted,rheap)

!<description>
  ! This routine attaches a directory-on-disc container to the external storage
  ! management.
!</description>

!<input>
  ! Path to the directory that saves the data.
  ! The directory must exist. If this is "", the current directory is used.
  character(LEN=*), intent(IN) :: spath

  ! OPTIONAL: If set to TRUE, the data in the container will be saved
  ! in a formatted file format. If set to FALSE (default), the data will
  ! be saved unformatted (which is machine dependent but faster).
  logical, intent(IN), optional :: bformatted
  
  ! OPTIONAL: Maximum size of the container (in Megabytes).
  ! -1 or not present = arbitrary.
  integer, intent(IN), optional :: imaxMem

  ! OPTIONAL: Basic filename of files that are stored on disc. The files will get
  ! the name "[filename].[handlenr]". If not specified, a default filename
  ! will be used.
  character(LEN=*), intent(IN), optional :: sfilename
  
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is initialised.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_exstorageContainer), dimension(:), pointer :: p_RstorageContainers
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Pointer to the heap to initialise
    type(t_exstorageBlock), pointer :: p_rheap

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    ! Create a new container.
    if (.not. associated(p_rheap%p_RstorageContainers)) then
      allocate(p_RstorageContainers(1))
    else
      allocate(p_RstorageContainers(size(p_rheap%p_RstorageContainers)+1))
      p_RstorageContainers(1:size(p_rheap%p_RstorageContainers)) = &
          p_rheap%p_RstorageContainers(:)
      
      deallocate(p_rheap%p_RstorageContainers)
    end if
    p_rheap%p_RstorageContainers => p_RstorageContainers
    
    p_rcontainer => p_rheap%p_RstorageContainers(size(p_rheap%p_RstorageContainers))
    
    ! Initialise the container as a Directory container
    p_rcontainer%ctype = EXSTOR_CONT_DIRECTORY
    if (spath .eq. '') then
      p_rcontainer%spath = './'
    else
      p_rcontainer%spath = trim(spath)//'/'
    end if
    if(present(imaxMem)) p_rcontainer%imaxStorageMB = imaxMem
    if(present(sfilename)) p_rcontainer%sfilename = sfilename
    if(present(bformatted)) p_rcontainer%bformatted = bformatted
    
  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_done(rheap)

!<description>
  ! This routine cleans up the storage management. All data on the
  ! heap is released from memory.
!</description>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is cleaned up.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap to initialise
    type(t_exstorageBlock), pointer :: p_rheap

    integer :: i,ihandle

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    ! Delete all data from the heap
    do i = 1,size(p_rheap%p_Rdescriptors)
      ! Don't pass i as handle as storage_free will set the handle
      ! passed to it to 0!
      ihandle = i
      if (p_rheap%p_Rdescriptors(i)%cdataType .ne. ST_NOHANDLE) &
        call exstor_free(ihandle,rheap)
    end do

    ! Clean up the memory management block
    p_rheap%nhandlesTotal = 0
    p_rheap%ihandlesDelta = 0
    p_rheap%p_inextFreeHandle = 0
    p_rheap%p_ilastFreeHandle = 0
    p_rheap%ihandlesInUse = 0

    ! Release the descriptors
    deallocate(p_rheap%p_IfreeHandles)
    deallocate(p_rheap%p_Rdescriptors)
    
    ! Release the storage containers
    deallocate(p_rheap%p_RstorageContainers)

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_info(bprintContainers, bprintHandles,rheap)

!<description>
  ! This routine prints information about the current memory consumption
  ! in a memory block to screen.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the information about the storage containers 
  ! is printed to the terminal.
  logical, intent(IN), optional :: bprintContainers

  ! OPTIONAL: If set to TRUE, the handles still remaining in the
  ! heap together with their names are printed to the terminal.
  logical, intent(IN), optional :: bprintHandles

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_exstorageBlock), intent(IN), target, optional :: rheap
!</input>

!</subroutine>

  ! local variables
  integer :: i

  ! Pointer to the heap
  type(t_exstorageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    call output_line ('External memory heap statistics:')
    call output_line ('--------------------------------')
    if (present(bprintHandles)) then
      if (bprintHandles .and. (p_rheap%ihandlesInUse .gt. 0)) then
        call output_line ('Handles on the heap: ')
        call output_lbrk ()
        ! Loop through the heap and search allocated handles
        do i=1,size(p_rheap%p_IfreeHandles)
          if (p_rheap%p_Rdescriptors(i)%cdataType .ne. ST_NOHANDLE) then
            if (p_rheap%p_Rdescriptors(i)%idimension .eq. 1) then
              call output_line ( &
                   'Handle ' // trim(sys_siL(i,10)) // ', 1D, Length=' // &
                   trim(sys_siL(int(p_rheap%p_Rdescriptors(i)%dmemBytes),15)) //&
                   ', Type=' // trim(sys_siL(p_rheap%p_Rdescriptors(i)%cdataType,15)) //&
                   ' Name=' // trim(adjustl(p_rheap%p_Rdescriptors(i)%sname)) )
            else
              call output_line ( &
                   'Handle ' // trim(sys_siL(i,10)) // ', 2D, Length=' // &
                   trim(sys_siL(int(p_rheap%p_Rdescriptors(i)%dmemBytes),15)) // &
                   ', Type=' // trim(sys_siL(p_rheap%p_Rdescriptors(i)%cdataType,15)) //&
                   ' Name=' // trim(adjustl(p_rheap%p_Rdescriptors(i)%sname)) )
            end if
          end if
        end do
        call output_lbrk ()
      end if
    end if

    if (present(bprintContainers)) then
      if (bprintContainers) then
        call output_line ('Storage containers: ')
        
        do i=1,size(p_rheap%p_RstorageContainers)
          call output_lbrk ()
          
          call output_line ('Storage container:  ' // trim(sys_siL(i,10)))
          
          select case (p_rheap%p_RstorageContainers(i)%ctype)
          case (EXSTOR_CONT_RAMDRIVE)
            call output_line ('Type:               RamDrive')
          case (EXSTOR_CONT_DIRECTORY)
            call output_line ('Type:               Directory ('//&
                trim(p_rheap%p_RstorageContainers(i)%spath)//')')
          end select
          
          if (p_rheap%p_RstorageContainers(i)%imaxStorageMB .eq. -1) then
            call output_line ('Max. memory(MB):    infinite')
          else
            call output_line ('Max. memory(MB):    '//&
              trim(sys_siL(p_rheap%p_RstorageContainers(i)%imaxStorageMB,10)))
          end if
          
          call output_line ('Current memory(MB): '//&
              trim(sys_siL(int(p_rheap%p_RstorageContainers(i)%dcurrentStorage&
                   /1000000.0_DP,I32),10)))

        end do
      
        call output_lbrk ()
      end if
    end if
                      
    call output_line ('Number of storage containers:    '//&
                      trim(sys_siL(size(p_rheap%p_RstorageContainers),15)))
    call output_line ('Number of handles in use:        '//&
                      trim(sys_siL(p_rheap%ihandlesInUse,15)))
    if (p_rheap%dtotalMem .gt. real(huge(0),DP)) then
      call output_line ('Memory in use (bytes):           '//&
                        trim(sys_sdL(p_rheap%dtotalMem,0)))
    else
      call output_line ('Memory in use (bytes):           '//&
                        trim(sys_siL(int(p_rheap%dtotalMem),15)))
    end if
    call output_line ('Current total number of handles: '//&
                      trim(sys_siL(size(p_rheap%p_IfreeHandles),15)))
    call output_line ('Maximum number of handles used:  '//&
                      trim(sys_siL(p_rheap%nhandlesInUseMax,15)))

    if (p_rheap%dtotalMem .gt. real(huge(0),DP)) then
      call output_line ('Maximum used memory (bytes):     '//&
                        trim(sys_sdL(p_rheap%dtotalMemMax,0)))
    else
      call output_line ('Maximum used memory (bytes):     '//&
                        trim(sys_siL(int(p_rheap%dtotalMemMax),15)))
    end if
  end subroutine

!************************************************************************

!<function>

  integer function exstor_newhandle (rheap) result(ihandle)

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
  type(t_exstorageBlock), intent(INOUT) :: rheap

!</inputoutput>

!</function>

    ! local variables
    type(t_exstorageNode), dimension(:), pointer :: p_Rdescriptors => null()
    integer, dimension(:), pointer :: p_IfreeHandles => null()
    integer :: i

    if (rheap%nhandlesTotal .le. 0) then
      call output_line ('Heap not initialised!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'exstor_newhandle')
      call sys_halt()
    end if

    ! Handles available?

    if (rheap%ihandlesInUse .ge. rheap%nhandlesTotal) then

      ! All handles are in use. We have to modify our ring to accept more
      ! handles.
      !
      ! At first, reallocate the descriptor-array and the queue-array with
      ! the new size.
      allocate (p_Rdescriptors (rheap%nhandlesTotal + rheap%ihandlesDelta) )
      allocate (p_IfreeHandles (rheap%nhandlesTotal + rheap%ihandlesDelta) )

      ! Copy the content, release the old arrays and replace them by the new
      ! ones.
      p_Rdescriptors(1:rheap%nhandlesTotal) = rheap%p_Rdescriptors(1:rheap%nhandlesTotal)
      p_IfreeHandles(1:rheap%nhandlesTotal) = rheap%p_IfreeHandles(1:rheap%nhandlesTotal)

      deallocate(rheap%p_Rdescriptors)
      deallocate(rheap%p_IfreeHandles)
      rheap%p_Rdescriptors => p_Rdescriptors
      rheap%p_IfreeHandles => p_IfreeHandles

      ! Add the new handles to the list of 'free' handles.
      do i=rheap%nhandlesTotal+1, rheap%nhandlesTotal + rheap%ihandlesDelta
        p_IfreeHandles (i) = i
      end do

      ! The first new 'free' handle is not at position...
      rheap%p_inextFreeHandle = rheap%nhandlesTotal+1

      ! And the last 'free' handle is at the end of the new list.
      rheap%p_ilastFreeHandle = rheap%nhandlesTotal + rheap%ihandlesDelta

      ! Modify the heap structure - we have more handles now.
      rheap%nhandlesTotal = rheap%nhandlesTotal + rheap%ihandlesDelta

    end if

    ! Get the new handle...
    ihandle = rheap%p_IfreeHandles (rheap%p_inextFreeHandle)

    ! and modify our queue pointers that we use a new one.
    rheap%p_inextFreeHandle = mod(rheap%p_inextFreeHandle,rheap%nhandlesTotal)+1

    rheap%ihandlesInUse = rheap%ihandlesInUse + 1

    rheap%nhandlesInUseMax = max(rheap%nhandlesInUseMax,rheap%ihandlesInUse)

  end function

!************************************************************************

!<subroutine>

  subroutine exstor_releasehandle (ihandle,rheap)

!<description>
  ! This routine releases a handle from the heap structure rheap.
  ! Memory is not deallocated, simply the structures are cleaned up.
!</description>

!<input>
  ! The handle to release
  integer, intent(INOUT) :: ihandle
!</input>

!<inputoutput>
  ! The heap structure where to release the handle from.
  type(t_exstorageBlock), intent(INOUT) :: rheap
!</inputoutput>

!</subroutine>

    type(t_exstorageNode), pointer :: p_rnode
    integer :: icontainer

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
    p_rnode%bcontainerBound = .false.
    p_rnode%cinitNewBlock = ST_NEWBLOCK_NOINIT
    p_rnode%Isize(1:2) = (/0,0/)

    ! Handle ihandle is available now - put it to the list of available handles.
    rheap%p_ilastFreeHandle = mod(rheap%p_ilastFreeHandle,rheap%nhandlesTotal) + 1
    rheap%p_IfreeHandles (rheap%p_ilastFreeHandle) = ihandle

    rheap%ihandlesInUse = rheap%ihandlesInUse - 1

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_free (ihandle, rheap)

!<description>
  ! This routine releases a handle from a heap and deallocates the
  ! associated memory. ihandle is set to ST_NOHANDLE upon return.
!</description>

!<inputoutput>

  ! Handle of the memory block to be releases
  integer :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!</subroutine>

  ! local variables

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Releasing ST_NOHANDLE is not allowed!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_free')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Is the node associated at all?
    if (p_rnode%cdataType .eq. ST_NOHANDLE) then
      call output_line ('Trying to release nonexistent handle: '//sys_siL(ihandle,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_free')
      call sys_halt()
    end if

    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    ! Release the memory assigned to that handle. How this is released depends
    ! on the container type.
    select case(p_rcontainer%ctype)

    case (EXSTOR_CONT_RAMDRIVE)
      ! The handle identifies another handle in the main memory storage subsystem.
      ! Note that the handle may be ST_NOHANDLE -- which is the case if no data
      ! is associated to this up to now.
      if (p_rnode%istorageHandle .ne. ST_NOHANDLE) &
        call storage_free (p_rnode%istorageHandle)

    case (EXSTOR_CONT_DIRECTORY)
      ! The handle identifies a file in a directory. Delete the file.
      call io_deleteFile (trim(p_rcontainer%spath)//p_rnode%sfilename)

    end select

    ! Release the handle itself.
    call exstor_releasehandle (ihandle,p_rheap)

    ! And finally reset the handle to ST_NOHANDLE.
    ihandle = ST_NOHANDLE

  end subroutine

!************************************************************************

!<function>

  integer function exstor_getContainer (rheap,dsize)

!<description>
  ! This routine determines a container that is large enough to hold
  ! a memory block of size dsize.
!</description>

!<input>
  ! Local heap structure where a container should be searched.
  type(t_exstorageBlock), intent(IN), target, optional :: rheap
  
  ! Number of bytes (encoded as double) which are to be stored.
  real(DP), intent(IN) :: dsize
!</input>

!<result>
  ! An id of a container that is large enough to hold the memory.
  ! If no container is found, Id 1 (the standard Id for the RamDrive
  ! container) is returned.
!</result>

!</function>

    integer :: icontainer
    
    ! Which strategy should we use for searching for a memory block?
    select case (rheap%ccontainerStrategy)
    case (EXSTOR_STRAT_LASTCONTAINER)
      ! Start at the end of the container list and proceed to the first container
      ! until we find one that can hold the data.
      do icontainer = size(rheap%p_RstorageContainers),1,-1
        if ((rheap%p_RstorageContainers(icontainer)%imaxStorageMB .eq. -1) .or. &
            ((rheap%p_RstorageContainers(icontainer)%dcurrentStorage + &
              dsize)/1000000.0_DP .lt. &
              real(rheap%p_RstorageContainers(icontainer)%imaxStorageMB,DP))) then
          exstor_getContainer = icontainer
          return
        end if
      end do
    
      ! Return the RamDrive container
      exstor_getContainer = 1
    end select

  end function

!************************************************************************

!<subroutine>

  subroutine exstor_new1D (scall, sname, isize, ctype, ihandle, &
                           cinitNewBlock, icontainerId, rheap)

!<description>
  ! This routine reserves a 1D memory block of desired size and type.
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(IN) :: scall

  !clear name of data field
  character(LEN=*), intent(IN) :: sname

  !requested storage size
  integer(I32), intent(IN) :: isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT,ST_NEWBLOCK_ORDERED)
  integer, intent(IN) :: cinitNewBlock
  
  ! OPTIONAL: Id of the storage container (1,2,3,...) that should maintain 
  ! the memory block. If not specified, a container will automatically be 
  ! chosen (depending on the allocation strategy ccontainerStrategy of the
  ! container).
  ! Container-Id=1 always identifies the standard RamDrive container that
  ! saves data into the main memory. Container-Id=0 automatically determines
  ! a container (just as if the parameter is not specified).
  ! Note that if a valid container id <> 0 is specified, the memory block will
  ! be bound to that container and not be automatically moved to another one.
  integer, intent(IN), optional :: icontainerId

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(OUT) :: ihandle

!</output>

!</subroutine>

    ! local variables
    integer :: icontainer
    logical :: bbound
    real(DP) :: dmemsize
    character(LEN=SYS_NAMELEN) :: snameBackup

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode

    if (isize .eq. 0) then
      call output_line ('isize=0', &
                        OU_CLASS_WARNING,OU_MODE_STD,'exstor_new1D')
      ihandle = ST_NOHANDLE
      return
    end if

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if
    
    ! Get the size of the memory block
    select case (ctype)
    case (ST_SINGLE)
      dmemsize = real(isize,DP)*real(ST_SINGLE2BYTES)
    case (ST_DOUBLE)
      dmemsize = real(isize,DP)*real(ST_DOUBLE2BYTES)
    case (ST_INT)
      dmemsize = real(isize,DP)*real(ST_INT2BYTES)
    case (ST_LOGICAL)
      dmemsize = real(isize,DP)*real(ST_LOGICAL2BYTES)
    case (ST_CHAR)
      dmemsize = real(isize,DP)*real(ST_CHAR2BYTES)
    case DEFAULT
      call output_line ('Unknown mem type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_new1D')
      call sys_halt()
    end select
    
    ! Determine the container id.
    icontainer = 0
    bbound = .false.
    if (present(icontainerId)) then
      icontainer = icontainerId
      
      ! The memory block is bound top that container if the
      ! container id is specified.
      bbound = icontainer .ne. 0
    end if
    
    if (icontainer .eq. 0) then
      ! We have to find a container that is large enough to hold the data.
      ! Automatically determine a container.
      icontainer = exstor_getContainer (p_rheap,dmemsize)
    end if
    
    if ((icontainer .lt. 1) .or. &
        (icontainer .gt. size(p_rheap%p_RstorageContainers))) then
      call output_line ('Invalid container id!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_new1D')
      call sys_halt()
    end if

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
    if (p_rheap%dtotalMem .gt. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem
      
    ! Notify the container about the new memory
    p_rheap%p_RstorageContainers(icontainer)%dcurrentStorage = &
        p_rheap%p_RstorageContainers(icontainer)%dcurrentStorage +  p_rnode%dmemBytes
        
    ! Some container-specific initialisation
    select case (p_rheap%p_RstorageContainers(icontainer)%ctype)
    case (EXSTOR_CONT_DIRECTORY)
      ! Create a filename based on the filemane template of the container
      p_rnode%sfilename = &
          trim(p_rheap%p_RstorageContainers(icontainer)%sfilename)//'.'//&
          trim(sys_siL(ihandle,10))
    end select

    ! Note: This routine will not immediately allocate memory!
    ! the actual memory allocation (or file creation) is done
    ! the first time, the memory is addressed!

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_new2D (scall, sname, Isize, ctype, ihandle, &
                           cinitNewBlock, icontainerId, rheap)

!<description>
  ! This routine reserves a 1D memory block of desired size and type.
!</description>

!<input>

  ! name of the calling routine
  character(LEN=*), intent(IN) :: scall

  ! clear name of data field
  character(LEN=*), intent(IN) :: sname

  ! requested storage size; DIMENSION(2)
  integer(I32), dimension(:), intent(IN) :: Isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT,ST_NEWBLOCK_ORDERED)
  integer, intent(IN) :: cinitNewBlock
  
  ! OPTIONAL: Id of the storage container (1,2,3,...) that should maintain 
  ! the memory block. If not specified, a container will automatically be 
  ! chosen (depending on the allocation strategy ccontainerStrategy of the
  ! container).
  ! Container-Id=1 always identifies the standard RamDrive container that
  ! saves data into the main memory. Container-Id=0 automatically determines
  ! a container (just as if the parameter is not specified).
  ! Note that if a valid container id <> 0 is specified, the memory block will
  ! be bound to that container and not be automatically moved to another one.
  integer, intent(IN), optional :: icontainerId

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!<output>
  ! Handle of the memory block.
  integer, intent(OUT) :: ihandle
!</output>

!</subroutine>

    ! local variables
    integer :: icontainer
    logical :: bbound
    real(DP) :: dmemsize
    character(LEN=SYS_NAMELEN) :: snameBackup

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode

    if ((Isize(1) .eq. 0) .or. (Isize(2) .eq. 0)) then
      call output_line ('Isize=0', &
                        OU_CLASS_WARNING,OU_MODE_STD,'exstor_new2D')
      ihandle = ST_NOHANDLE
      return
    end if

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if
    
    ! Get the size of the memory block
    select case (ctype)
    case (ST_SINGLE)
      dmemsize = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_SINGLE2BYTES)
    case (ST_DOUBLE)
      dmemsize = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_DOUBLE2BYTES)
    case (ST_INT)
      dmemsize = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_INT2BYTES)
    case (ST_LOGICAL)
      dmemsize = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_LOGICAL2BYTES)
    case (ST_CHAR)
      dmemsize = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_CHAR2BYTES)
    case DEFAULT
      call output_line ('Unknown mem type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_new2D')
      call sys_halt()
    end select
    
    ! Determine the container id.
    icontainer = 0
    bbound = .false.
    if (present(icontainerId)) then
      icontainer = icontainerId
      
      ! The memory block is bound top that container if the
      ! container id is specified.
      bbound = icontainer .ne. 0
    end if
    
    if (icontainer .eq. 0) then
      ! We have to find a container that is large enough to hold the data.
      ! Automatically determine a container.
      icontainer = exstor_getContainer (p_rheap,dmemsize)
    end if
    
    if ((icontainerId .lt. 1) .or. &
        (icontainerId .gt. size(p_rheap%p_RstorageContainers))) then
      call output_line ('Invalid container id!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_new2D')
      call sys_halt()
    end if

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
    if (p_rheap%dtotalMem .gt. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem
      
    ! Notify the container about the new memory
    p_rheap%p_RstorageContainers(icontainer)%dcurrentStorage = &
        p_rheap%p_RstorageContainers(icontainer)%dcurrentStorage +  p_rnode%dmemBytes

    ! Some container-specific initialisation
    select case (p_rheap%p_RstorageContainers(icontainer)%ctype)
    case (EXSTOR_CONT_DIRECTORY)
      ! Create a filename based on the filemane template of the container
      p_rnode%sfilename = &
          trim(p_rheap%p_RstorageContainers(icontainer)%sfilename)//'.'//&
          trim(sys_siL(ihandle,10))
    end select

    ! Note: This routine will not immediately allocate memory!
    ! the actual memory allocation (or file creation) is done
    ! the first time, the memory is addressed!

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_getsize1D (ihandle, isize, rheap)

!<description>
  ! Returns the length of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block 
  integer, intent(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!<output>
  ! Length of the array identified by ihandle.
  integer(I32), intent(OUT) :: isize
!</output>

!</subroutine>

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getsize1D')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Is the node associated at all?
    if (p_rnode%cdataType .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid: '//sys_siL(ihandle,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getsize1D')
      call sys_halt()
    end if

    ! What are we?
    if (p_rnode%idimension .ne. 1) then
      call output_line ('Handle '//trim(sys_siL(ihandle,10))//' is not 1-dimensional!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getsize1D')
      call sys_halt()
    end if
    
    isize = p_rnode%Isize(1)

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_getsize2D (ihandle, isize, rheap)

!<description>
  ! Returns the length of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!<output>
  ! Length of each dimension of the array identified by ihandle.
  integer(I32), dimension(:), intent(OUT) :: Isize
!</output>

!</subroutine>

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Is the node associated at all?
    if (p_rnode%cdataType .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid: '//sys_siL(ihandle,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      call sys_halt()
    end if

    ! What are we?
    if (p_rnode%idimension .ne. 2) then
      call output_line ('Handle '//trim(sys_siL(ihandle,10))//' is not 2-dimensional!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      call sys_halt()
    end if
    
    Isize(1:2) = p_rnode%Isize(1:2)
    
  end subroutine
  
!************************************************************************

!<subroutine>

  subroutine exstor_getdata_int (ihandle, Iarray, rheap)

!<description>
  ! This routine reads data from an external storage container and
  ! saves it to a local array.
!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! An array that is filled with the data identified by the handle.
  integer(I32), dimension(:), intent(OUT) :: Iarray

!</output>

!</subroutine>

    ! local variables
    integer :: iorder

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%cdataType .ne. ST_INT) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    if (p_rnode%Isize(1) .ne. size(Iarray)) then
      call output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    ! If the memory block is a pre-initialised block, we can directly
    ! fill it with data...    
    select case (p_rnode%cinitNewBlock)
    case (ST_NEWBLOCK_ZERO)
      ! Fill the destination array with zero, that's all
      call lalg_clearVectorInt (Iarray)

    case (ST_NEWBLOCK_ORDERED)
      ! Fill the destination array by increasing numbers
      do iorder=1,size(Iarray)
        Iarray(iorder) = int(iorder,I32)
      end do

    case (ST_NEWBLOCK_NOINIT)
      ! Ok, there should be data behind. How to handle the data depends
      ! on the container type
      select case(p_rcontainer%ctype)
      
      case (EXSTOR_CONT_RAMDRIVE)
        ! This is a RamDrive container. We use the memory management of the
        ! storage.f90 to maintain it.
        call getdata_ramdrive (p_rnode,Iarray)
                
      case (EXSTOR_CONT_DIRECTORY)
        ! This is a directory container maintaining the data as files on the 
        ! hard disc.
        call getdata_directory (p_rcontainer,p_rnode,Iarray)

      end select
        
    end select

  contains
  
    ! -------------------------------------------------------------------------
    subroutine getdata_ramdrive (rnode,dataarray)
    
    ! Retrieves the data from a RamDrive container
    
    ! The storage block containing the data
    type(t_exstorageNode), intent(IN) :: rnode
    
    ! The destination array for the data
    integer(I32), dimension(:), intent(OUT) :: dataarray
    
      ! local variables
      integer(I32), dimension(:), pointer :: p_data

      ! At first: Do we have data at all?
      if (rnode%istorageHandle .eq. ST_NOHANDLE) then
        call output_line ('Trying to read noninitialised memory! Handle: '//&
            sys_siL(ihandle,10), OU_CLASS_ERROR,OU_MODE_STD,'getdata_ramdrive')
        call sys_halt()
      end if
    
      ! Get the memory block and copy it.
      call storage_getbase_int (rnode%istoragehandle,p_data)
      call lalg_copyVectorInt (p_data,dataarray)
      
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine getdata_directory (rcontainer,rnode,dataarray)
    
    ! Retrieves the data from a directory container
    
    ! The storage container assigned to the storage block
    type(t_exstorageContainer), intent(IN) :: rcontainer
        
    ! The storage block containing the data
    type(t_exstorageNode), intent(IN) :: rnode
    
    ! The destination array for the data
    integer(I32), dimension(:), intent(OUT) :: dataarray
    
      ! local variables
      integer :: cf
      
      ! Open the file to read data
      call io_openFileForReading(trim(rcontainer%spath)//rnode%sfilename, &
          cf, bformatted=rcontainer%bformatted)
      
      ! Read the data from the file
      if (rcontainer%bformatted) then
        read(cf,*) dataarray(:)
      else
        read(cf) dataarray(:)
      end if
      
      ! Close the file, finish
      close (cf)
      
    end subroutine

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_setdata_int (ihandle, Iarray, rheap)

!<description>
  ! This routine writes data from a local array top an external storage 
  ! container.
!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

  ! An array with data that should be stored on the external storage
  ! container.
  integer(I32), dimension(:), intent(IN) :: Iarray
!</input>

!</subroutine>

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%cdataType .ne. ST_INT) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    if (p_rnode%Isize(1) .ne. size(Iarray)) then
      call output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    ! The data block will now contain data.
    p_rnode%cinitNewBlock = ST_NEWBLOCK_NOINIT

    ! How to handle the data depends on the container type.
    select case(p_rcontainer%ctype)
    
    case (EXSTOR_CONT_RAMDRIVE)
      ! This is a RamDrive container. We use the memory management of the
      ! storage.f90 to maintain it.
      call setdata_ramdrive (p_rnode,Iarray)
              
    case (EXSTOR_CONT_DIRECTORY)
      ! This is a directory container maintaining the data as files on the 
      ! hard disc.
      call setdata_directory (p_rcontainer,p_rnode,Iarray)

    end select
        
  contains
  
    ! -------------------------------------------------------------------------
    subroutine setdata_ramdrive (rnode,dataarray)
    
    ! Writes data to a RamDrive container
    
    ! The storage block containing the data
    type(t_exstorageNode), intent(INOUT) :: rnode
    
    ! The source array with the data
    integer(I32), dimension(:), intent(IN) :: dataarray
    
      ! local variables
      integer(I32), dimension(:), pointer :: p_data

      ! Copy data with the storage-copy command. If necessary, new memory
      ! is allocated.
      if (rnode%istorageHandle .eq. ST_NOHANDLE) then
        ! We have to allocate memory.
        call storage_new ('setdata_ramdrive','exdata',rnode%Isize(1),&
            rnode%cdataType,rnode%istorageHandle,ST_NEWBLOCK_NOINIT)
      end if
    
      ! Get the memory block and copy it.
      call storage_getbase_int (rnode%istoragehandle,p_data)
      call lalg_copyVectorInt (dataarray,p_data)
      
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine setdata_directory (rcontainer,rnode,dataarray)
    
    ! Writes data to a directory container
    
    ! The storage container assigned to the storage block
    type(t_exstorageContainer), intent(IN) :: rcontainer
        
    ! The storage block containing the data
    type(t_exstorageNode), intent(INOUT) :: rnode
    
    ! The source array with the data
    integer(I32), dimension(:), intent(IN) :: dataarray
    
      ! local variables
      integer :: cf
      
      ! Open the file to write data
      call io_openFileForWriting(trim(rcontainer%spath)//rnode%sfilename, &
          cf, SYS_REPLACE, bformatted=rcontainer%bformatted)
      
      ! Write the data to the file
      if (rcontainer%bformatted) then
        write(cf,*) dataarray(:)
      else
        write(cf) dataarray(:) 
      end if
      
      ! Close the file, finish
      close (cf)
      
    end subroutine

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_getdata_single (ihandle, Farray, rheap)

!<description>
  ! This routine reads data from an external storage container and
  ! saves it to a local array.
!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! An array that is filled with the data identified by the handle.
  real(SP), dimension(:), intent(OUT) :: Farray

!</output>

!</subroutine>

    ! local variables
    integer :: iorder

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%cdataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    if (p_rnode%Isize(1) .ne. size(Farray)) then
      call output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    ! If the memory block is a pre-initialised block, we can directly
    ! fill it with data...    
    select case (p_rnode%cinitNewBlock)
    case (ST_NEWBLOCK_ZERO)
      ! Fill the destination array with zero, that's all
      call lalg_clearVectorSngl (Farray)

    case (ST_NEWBLOCK_ORDERED)
      ! Fill the destination array by increasing numbers
      do iorder=1,size(Farray)
        Farray(iorder) = real(iorder,SP)
      end do

    case (ST_NEWBLOCK_NOINIT)
      ! Ok, there should be data behind. How to handle the data depends
      ! on the container type
      select case(p_rcontainer%ctype)
      
      case (EXSTOR_CONT_RAMDRIVE)
        ! This is a RamDrive container. We use the memory management of the
        ! storage.f90 to maintain it.
        call getdata_ramdrive (p_rnode,Farray)
                
      case (EXSTOR_CONT_DIRECTORY)
        ! This is a directory container maintaining the data as files on the 
        ! hard disc.
        call getdata_directory (p_rcontainer,p_rnode,Farray)

      end select
        
    end select

  contains
  
    ! -------------------------------------------------------------------------
    subroutine getdata_ramdrive (rnode,dataarray)
    
    ! Retrieves the data from a RamDrive container
    
    ! The storage block containing the data
    type(t_exstorageNode), intent(IN) :: rnode
    
    ! The destination array for the data
    real(SP), dimension(:), intent(OUT) :: dataarray
    
      ! local variables
      real(SP), dimension(:), pointer :: p_data

      ! At first: Do we have data at all?
      if (rnode%istorageHandle .eq. ST_NOHANDLE) then
        call output_line ('Trying to read noninitialised memory! Handle: '//&
            sys_siL(ihandle,10), OU_CLASS_ERROR,OU_MODE_STD,'getdata_ramdrive')
        call sys_halt()
      end if
    
      ! Get the memory block and copy it.
      call storage_getbase_single (rnode%istoragehandle,p_data)
      call lalg_copyVectorSngl (p_data,dataarray)
      
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine getdata_directory (rcontainer,rnode,dataarray)
    
    ! Retrieves the data from a directory container
    
    ! The storage container assigned to the storage block
    type(t_exstorageContainer), intent(IN) :: rcontainer
        
    ! The storage block containing the data
    type(t_exstorageNode), intent(IN) :: rnode
    
    ! The destination array for the data
    real(SP), dimension(:), intent(OUT) :: dataarray
    
      ! local variables
      integer :: cf
      
      ! Open the file to read data
      call io_openFileForReading(trim(rcontainer%spath)//rnode%sfilename, &
          cf, bformatted=rcontainer%bformatted)
      
      ! Read the data from the file
      if (rcontainer%bformatted) then
        read(cf,*) dataarray(:)
      else
        read(cf) dataarray(:)
      end if
      
      ! Close the file, finish
      close (cf)
      
    end subroutine

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_setdata_single (ihandle, Farray, rheap)

!<description>
  ! This routine writes data from a local array top an external storage 
  ! container.
!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

  ! An array with data that should be stored on the external storage
  ! container.
  real(SP), dimension(:), intent(IN) :: Farray
!</input>

!</subroutine>

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%cdataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    if (p_rnode%Isize(1) .ne. size(Farray)) then
      call output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    ! The data block will now contain data.
    p_rnode%cinitNewBlock = ST_NEWBLOCK_NOINIT

    ! How to handle the data depends on the container type.
    select case(p_rcontainer%ctype)
    
    case (EXSTOR_CONT_RAMDRIVE)
      ! This is a RamDrive container. We use the memory management of the
      ! storage.f90 to maintain it.
      call setdata_ramdrive (p_rnode,Farray)
              
    case (EXSTOR_CONT_DIRECTORY)
      ! This is a directory container maintaining the data as files on the 
      ! hard disc.
      call setdata_directory (p_rcontainer,p_rnode,Farray)

    end select
        
  contains
  
    ! -------------------------------------------------------------------------
    subroutine setdata_ramdrive (rnode,dataarray)
    
    ! Writes data to a RamDrive container
    
    ! The storage block containing the data
    type(t_exstorageNode), intent(INOUT) :: rnode
    
    ! The source array with the data
    real(SP), dimension(:), intent(IN) :: dataarray
    
      ! local variables
      real(SP), dimension(:), pointer :: p_data

      ! Copy data with the storage-copy command. If necessary, new memory
      ! is allocated.
      if (rnode%istorageHandle .eq. ST_NOHANDLE) then
        ! We have to allocate memory.
        call storage_new ('setdata_ramdrive','exdata',rnode%Isize(1),&
            rnode%cdataType,rnode%istorageHandle,ST_NEWBLOCK_NOINIT)
      end if
    
      ! Get the memory block and copy it.
      call storage_getbase_single (rnode%istoragehandle,p_data)
      call lalg_copyVectorSngl (dataarray,p_data)
      
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine setdata_directory (rcontainer,rnode,dataarray)
    
    ! Writes data to a directory container
    
    ! The storage container assigned to the storage block
    type(t_exstorageContainer), intent(IN) :: rcontainer
        
    ! The storage block containing the data
    type(t_exstorageNode), intent(INOUT) :: rnode
    
    ! The source array with the data
    real(SP), dimension(:), intent(IN) :: dataarray
    
      ! local variables
      integer :: cf
      
      ! Open the file to write data
      call io_openFileForWriting(trim(rcontainer%spath)//rnode%sfilename, &
          cf, SYS_REPLACE, bformatted=rcontainer%bformatted)
      
      ! Write the data to the file
      if (rcontainer%bformatted) then
        write(cf,*) dataarray(:)
      else
        write(cf) dataarray(:) 
      end if
      
      ! Close the file, finish
      close (cf)
      
    end subroutine

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_getdata_double (ihandle, Darray, rheap)

!<description>
  ! This routine reads data from an external storage container and
  ! saves it to a local array.
!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! An array that is filled with the data identified by the handle.
  real(DP), dimension(:), intent(OUT) :: Darray

!</output>

!</subroutine>

    ! local variables
    integer :: iorder

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%cdataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    if (p_rnode%Isize(1) .ne. size(Darray)) then
      call output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_int')
      call sys_halt()
    end if

    ! If the memory block is a pre-initialised block, we can directly
    ! fill it with data...    
    select case (p_rnode%cinitNewBlock)
    case (ST_NEWBLOCK_ZERO)
      ! Fill the destination array with zero, that's all
      call lalg_clearVectorDble (Darray)

    case (ST_NEWBLOCK_ORDERED)
      ! Fill the destination array by increasing numbers
      do iorder=1,size(Darray)
        Darray(iorder) = real(iorder,DP)
      end do

    case (ST_NEWBLOCK_NOINIT)
      ! Ok, there should be data behind. How to handle the data depends
      ! on the container type
      select case(p_rcontainer%ctype)
      
      case (EXSTOR_CONT_RAMDRIVE)
        ! This is a RamDrive container. We use the memory management of the
        ! storage.f90 to maintain it.
        call getdata_ramdrive (p_rnode,Darray)
                
      case (EXSTOR_CONT_DIRECTORY)
        ! This is a directory container maintaining the data as files on the 
        ! hard disc.
        call getdata_directory (p_rcontainer,p_rnode,Darray)

      end select
        
    end select

  contains
  
    ! -------------------------------------------------------------------------
    subroutine getdata_ramdrive (rnode,dataarray)
    
    ! Retrieves the data from a RamDrive container
    
    ! The storage block containing the data
    type(t_exstorageNode), intent(IN) :: rnode
    
    ! The destination array for the data
    real(DP), dimension(:), intent(OUT) :: dataarray
    
      ! local variables
      real(DP), dimension(:), pointer :: p_data

      ! At first: Do we have data at all?
      if (rnode%istorageHandle .eq. ST_NOHANDLE) then
        call output_line ('Trying to read noninitialised memory! Handle: '//&
            sys_siL(ihandle,10), OU_CLASS_ERROR,OU_MODE_STD,'getdata_ramdrive')
        call sys_halt()
      end if
    
      ! Get the memory block and copy it.
      call storage_getbase_double (rnode%istoragehandle,p_data)
      call lalg_copyVectorDble (p_data,dataarray)
      
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine getdata_directory (rcontainer,rnode,dataarray)
    
    ! Retrieves the data from a directory container
    
    ! The storage container assigned to the storage block
    type(t_exstorageContainer), intent(IN) :: rcontainer

    ! The storage block containing the data
    type(t_exstorageNode), intent(IN) :: rnode
    
    ! The destination array for the data
    real(DP), dimension(:), intent(OUT) :: dataarray
    
      ! local variables
      integer :: cf
      
      ! Open the file to read data
      call io_openFileForReading(trim(rcontainer%spath)//rnode%sfilename, &
          cf, bformatted=rcontainer%bformatted)
      
      ! Read the data from the file
      if (rcontainer%bformatted) then
        read(cf,*) dataarray(:)
      else
        read(cf) dataarray(:)
      end if
      
      ! Close the file, finish
      close (cf)
      
    end subroutine

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_setdata_double (ihandle, Darray, rheap)

!<description>
  ! This routine writes data from a local array top an external storage 
  ! container.
!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

  ! An array with data that should be stored on the external storage
  ! container.
  real(DP), dimension(:), intent(IN) :: Darray
!</input>

!</subroutine>

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%cdataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    if (p_rnode%Isize(1) .ne. size(Darray)) then
      call output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_int')
      call sys_halt()
    end if

    ! The data block will now contain data.
    p_rnode%cinitNewBlock = ST_NEWBLOCK_NOINIT

    ! How to handle the data depends on the container type.
    select case(p_rcontainer%ctype)
    
    case (EXSTOR_CONT_RAMDRIVE)
      ! This is a RamDrive container. We use the memory management of the
      ! storage.f90 to maintain it.
      call setdata_ramdrive (p_rnode,Darray)
              
    case (EXSTOR_CONT_DIRECTORY)
      ! This is a directory container maintaining the data as files on the 
      ! hard disc.
      call setdata_directory (p_rcontainer,p_rnode,Darray)

    end select
        
  contains
  
    ! -------------------------------------------------------------------------
    subroutine setdata_ramdrive (rnode,dataarray)
    
    ! Writes data to a RamDrive container
    
    ! The storage block containing the data
    type(t_exstorageNode), intent(INOUT) :: rnode
    
    ! The source array with the data
    real(DP), dimension(:), intent(IN) :: dataarray
    
      ! local variables
      real(DP), dimension(:), pointer :: p_data

      ! Copy data with the storage-copy command. If necessary, new memory
      ! is allocated.
      if (rnode%istorageHandle .eq. ST_NOHANDLE) then
        ! We have to allocate memory.
        call storage_new ('setdata_ramdrive','exdata',rnode%Isize(1),&
            rnode%cdataType,rnode%istorageHandle,ST_NEWBLOCK_NOINIT)
      end if
    
      ! Get the memory block and copy it.
      call storage_getbase_double (rnode%istoragehandle,p_data)
      call lalg_copyVectorDble (dataarray,p_data)
      
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine setdata_directory (rcontainer,rnode,dataarray)
    
    ! Writes data to a directory container
    
    ! The storage container assigned to the storage block
    type(t_exstorageContainer), intent(IN) :: rcontainer
        
    ! The storage block containing the data
    type(t_exstorageNode), intent(INOUT) :: rnode
    
    ! The source array with the data
    real(DP), dimension(:), intent(IN) :: dataarray
    
      ! local variables
      integer :: cf
      
      ! Open the file to write data
      call io_openFileForWriting(trim(rcontainer%spath)//rnode%sfilename, &
          cf, SYS_REPLACE, bformatted=rcontainer%bformatted)
      
      ! Write the data to the file
      if (rcontainer%bformatted) then
        write(cf,*) dataarray(:)
      else
        write(cf) dataarray(:) 
      end if
      
      ! Close the file, finish
      close (cf)
      
    end subroutine

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_getdata_storage (ihandle, istoragehandle, rheap, rstorageheap)

!<description>
  ! This routine reads data from an external storage container and
  ! saves it to a local array in the storage management of the main
  ! memory.
!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle
  
  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle to a memory block in the main memory, maintained by the 
  ! storage.f90. Data read from the data container is directly saved
  ! to this memory block.
  ! If this is ST_NOHANDLE, a new handle is automatically created
  ! and the data is saved to it.
  integer, intent(INOUT) :: istoragehandle

  ! OPTIONAL: local heap structure of the storage management for the
  ! main memory. If not given, the global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rstorageheap
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ctype
    integer(I32) :: isize

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    integer(I32), dimension(:), pointer :: p_Idata

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
      call sys_halt()
    end if

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    if (istoragehandle .eq. ST_NOHANDLE) then
    
      ! Allocate a new memory block with the correct data type that is just
      ! as large that it can hold our data.
      
      select case(p_rnode%idimension)
      case (1)
        call storage_new ('exstor_getdata_storage','datacopy',p_rnode%Isize(1),&
            p_rnode%cdataType,istoragehandle,ST_NEWBLOCK_NOINIT,rstorageheap)
      case (2)
        call storage_new ('exstor_getdata_storage','datacopy',p_rnode%Isize,&
            p_rnode%cdataType,istoragehandle,ST_NEWBLOCK_NOINIT,rstorageheap)
      end select
      
    else
      call storage_getdatatype(istoragehandle,ctype)
      if (p_rnode%cdataType .ne. ctype) then
        call output_line ('Wrong data format!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
        call sys_halt()
      end if
      
      call storage_getsize (istoragehandle,isize)
      if (p_rnode%Isize(1) .ne. isize) then
        call output_line ('Data array has the wrong size!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
        call sys_halt()
      end if

    end if
    
    ! If the memory block is a pre-initialised block, we can directly
    ! fill it with data...    
    if (p_rnode%cinitNewBlock .ne. ST_NEWBLOCK_NOINIT) then
      call storage_initialiseBlock (istoragehandle, p_rnode%cinitNewBlock)
      return
    end if
    
    ! Ok, we have to copy data from the external storage to the memory block.
    ! This now depends on the data type...
    select case (p_rnode%cdataType)
    case (ST_DOUBLE)
      call storage_getbase_double (istoragehandle,p_Ddata,rstorageheap)
      call exstor_getdata_double (ihandle,p_Ddata)
    case (ST_SINGLE)
      call storage_getbase_single (istoragehandle,p_Fdata,rstorageheap)
      call exstor_getdata_single (ihandle,p_Fdata)
    case (ST_INT)
      call storage_getbase_int (istoragehandle,p_Idata,rstorageheap)
      call exstor_getdata_int (ihandle,p_Idata)
    case DEFAULT
      call output_line ('Unsupported data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
    end select

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_setdata_storage (ihandle, istoragehandle, rheap, rstorageheap)

!<description>
  ! This routine writes data from a local array in the storage management
  ! of the main memory to an external storage container.
!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle
  
  ! Handle to a memory block in the main memory, maintained by the 
  ! storage.f90. The data in this memory block is directly saved
  ! to the the data container.
  integer, intent(IN) :: istoragehandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap

  ! OPTIONAL: local heap structure of the storage management for the
  ! main memory. If not given, the global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rstorageheap
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ctype
    integer(I32) :: isize

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    type(t_exstorageContainer), pointer :: p_rcontainer
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    integer(I32), dimension(:), pointer :: p_Idata

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
      call sys_halt()
    end if

    if (istoragehandle .eq. ST_NOHANDLE) then
      call output_line ('istoragehandle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
      call sys_halt()
    end if

    ! Get the container and the storage node
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    p_rcontainer => p_rheap%p_RstorageContainers(p_rnode%icontainerId)
    
    call storage_getdatatype(istoragehandle,ctype)
    if (p_rnode%cdataType .ne. ctype) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
      call sys_halt()
    end if
    
    call storage_getsize (istoragehandle,isize)
    if (p_rnode%Isize(1) .ne. isize) then
      call output_line ('Data array has the wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
      call sys_halt()
    end if
    
    ! Ok, we have to copy data from the external storage to the memory block.
    ! This now depends on the data type...
    select case (p_rnode%cdataType)
    case (ST_DOUBLE)
      call storage_getbase_double (istoragehandle,p_Ddata)
      call exstor_setdata_double (ihandle,p_Ddata)
    case (ST_SINGLE)
      call storage_getbase_single (istoragehandle,p_Fdata)
      call exstor_setdata_single (ihandle,p_Fdata)
    case (ST_INT)
      call storage_getbase_int (istoragehandle,p_Idata)
      call exstor_setdata_int (ihandle,p_Idata)
    case DEFAULT
      call output_line ('Unsupported data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_setdata_storage')
    end select

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_clear (ihandle, rheap)

!<description>
  ! This routine clears an array identified by ihandle; all entries are
  ! overwritten by 0.
!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure. If not given, the global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!</subroutine>

    ! Pointer to the heap
    type(t_exstorageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbaseexternal
    end if

    if (ihandle .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exstor_clear')
      call sys_halt()
    end if

    ! Set the initialised-flag to ST_NEWBLOCK_ZERO. The data is
    ! trated as zero when it's accessed the next time.
    p_rheap%p_Rdescriptors(ihandle)%cinitNewBlock = ST_NEWBLOCK_ZERO

  end subroutine

!************************************************************************

!<subroutine>

  subroutine exstor_copy(h_source, h_dest, icontainerId, rheap)

!<description>
  ! This routine copies the information of one array to another.
  ! The structure of the arrays behind h_source and h_dest must be the
  ! same!
!</description>

!<input>
  ! Handle of the source array to copy
  integer, intent(IN) :: h_source

  ! OPTIONAL: If h_dest=ST_NOHANDLE, this parameter allows to specify
  ! a container id where new data is stored to.
  ! If h_dest<>ST_NOHANDLE, this parameter is ignored.
  integer, intent(IN), optional :: icontainerId

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_exstorageBlock), intent(INOUT), target, optional :: rheap
  
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it. Such a new
  ! block is allocated in the container icontainerId
  integer, intent(INOUT) :: h_dest
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_exstorageBlock), pointer :: p_rheap
    type(t_exstorageNode), pointer :: p_rnode
    integer :: ihandle
    
    ! Read in the data, create a temp array in the main memory.
    ihandle = ST_NOHANDLE
    call exstor_getdata_storage (h_source, ihandle, rheap)
    
    ! If necessary, allocate a new block.
    if (h_dest .eq. ST_NOHANDLE) then

      if(present(rheap)) then
        p_rheap => rheap
      else
        p_rheap => rbaseexternal
      end if

      if (h_source .eq. ST_NOHANDLE) then
        call output_line ('Handle invalid!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'exstor_getdata_storage')
        call sys_halt()
      end if
    
      ! Get the container and the storage node
      p_rnode => p_rheap%p_Rdescriptors(ihandle)
      
      ! Allocate a new memory block with the correct data type that is just
      ! as large that it can hold our data.
      
      select case(p_rnode%idimension)
      case (1)
        call exstor_new ('exstor_copy','datacopy',p_rnode%Isize(1),&
            p_rnode%cdataType,h_dest,ST_NEWBLOCK_NOINIT,icontainerId,rheap)
      case (2)
        call exstor_new ('exstor_copy','datacopy',p_rnode%Isize,&
            p_rnode%cdataType,h_dest,ST_NEWBLOCK_NOINIT,icontainerId,rheap)
      end select
      
    end if
    
    ! Write the data to a new block
    call exstor_setdata_storage (h_dest, ihandle, rheap)
    
    ! Release the temp memory
    call storage_free (ihandle)

  end subroutine

end module
