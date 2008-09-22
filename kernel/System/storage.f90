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
!# With 'storage_realloc', memory can be reallocated, that is, the size
!# of the memory block associated with a given handle is modified.
!# With 'storage_getbase_xxxx', the pointer corresponding to a handle can be
!# obtained.
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
!# 13.) storage_getname
!#      -> Get the name of an array on the heap.
!#
!# 14.) storage_realloc
!#      -> Reallocate a 1D or 2D array (only 2nd dimension)
!#
!# 15.) storage_initialiseBlock
!#      -> Initialise a storage block with zero (like storage_clear)
!#         or increasing number.
!#
!# 16.) storage_isEqual
!#      -> Checks if the content of two different handles is equal
!#
!# 17.) storage_isFree
!#      -> Checks if the given handle is free
!#
!# 18.) storage_getBlock
!#      -> Get the internal data of the storage block
!# </purpose>
!##############################################################################

module storage
  
  use fpersistaux
  use fsystem
  use genoutput
  use linearalgebra
  use uuid

  implicit none

!<constants>

  !<constantblock description="Storage block type identifiers">

  ! defines a non-allocated storage block handle
  integer, parameter :: ST_NOHANDLE = 0

  ! storage block contains single floats
  integer, parameter :: ST_SINGLE = 1

  ! storage block contains double floats
  integer, parameter :: ST_DOUBLE = 2

  ! storage block contains ints
  integer, parameter :: ST_INT = 3

  ! storage block contains logicals
  integer, parameter :: ST_LOGICAL = 4

  ! storage block contains characters
  integer, parameter :: ST_CHAR = 5

  !</constantblock>

  !<constantblock description="Constants for initialisation of memory on allocation">

  ! init new storage block with zeros (or .FALSE. for logicals)
  integer, parameter :: ST_NEWBLOCK_ZERO = 0

  ! no init new storage block
  integer, parameter :: ST_NEWBLOCK_NOINIT = 1

  ! init new storage block with 1,2,3,...,n
  ! Note: This initialization constant must NOT be used for initialisation of logicals
  !       or characters!!!
  integer, parameter :: ST_NEWBLOCK_ORDERED = 2

  !</constantblock>

  !<constantblock description="Constants for calculating memory">

  ! How many bytes has an integer?
  integer :: ST_INT2BYTES = I32

  ! How many bytes has a standard real?
  integer :: ST_SINGLE2BYTES = SP

  ! How many bytes has a double precision real?
  integer :: ST_DOUBLE2BYTES = DP

  ! How many bytes has a logical?
  integer :: ST_LOGICAL2BYTES = I32

  ! How many bytes has a character?
  ! Note: We are not 100% sure, but this may differ on other architectures... O_o
  integer :: ST_CHAR2BYTES = 1

  !</constantblock>

!</constants>

!<types>

  !<typeblock>

  ! Type block for describing a handle. This collects the number of the
  ! handle, the storage amount associated to it, the pointer to the memory
  ! location etc.

  type t_storageNode

    private

    ! Type of data associated to the handle (ST_NOHANDLE, ST_SINGLE,
    ! ST_DOUBLE, ST_INT)
    integer :: idataType = ST_NOHANDLE

    ! Dimension associated to the handle (0=not assigned, 1=1D, 2=2D array)
    integer :: idimension = 0

    ! The name of the array that is associated to that handle
    character(LEN=SYS_NAMELEN) :: sname = ''

    ! The name of the calling routine that modified this handle
    character(LEN=SYS_NAMELEN) :: scall = ''

    ! Amount of memory (in bytes) associated to this block.
    ! We store that as a double to allow storing numbers > 2GB !
    real(DP) :: dmemBytes = 0.0_DP

    ! Pointer to 1D real array or NULL() if not assigned
    real(SP), dimension(:), pointer       :: p_Fsingle1D    => null()

    ! Pointer to 1D double precision array or NULL() if not assigned
    real(DP), dimension(:), pointer       :: p_Ddouble1D    => null()

    ! Pointer to 1D integer array or NULL() if not assigned
    integer(I32), dimension(:), pointer   :: p_Iinteger1D   => null()

    ! Pointer to 1D logical array or NULL() if not assigned
    logical, dimension(:), pointer        :: p_Blogical1D   => null()

    ! Pointer to 1D character array or NULL() if not assigned
    character, dimension(:), pointer      :: p_Schar1D      => null()

    ! Pointer to 2D real array or NULL() if not assigned
    real(SP), dimension(:,:), pointer     :: p_Fsingle2D    => null()

    ! Pointer to 2D double precision array or NULL() if not assigned
    real(DP), dimension(:,:), pointer     :: p_Ddouble2D    => null()

    ! Pointer to 2D integer array or NULL() if not assigned
    integer(I32), dimension(:,:), pointer :: p_Iinteger2D   => null()

    ! Pointer to 2D logical array or NULL() if not assigned
    logical, dimension(:,:), pointer        :: p_Blogical2D => null()

    ! Pointer to 2D character array or NULL() if not assigned
    character, dimension(:,:), pointer      :: p_Schar2D    => null()

  end type t_storageNode

  !</typeblock>

  !<typeblock>

  ! This block represents a heap that maintains single, double precision
  ! and integer data. It contains a list of t_storageNode elements for all
  ! the handles.
  ! There's one global object of this type for the global storage management,
  ! but if necessary, an algorithm can create such a block locally, too,
  ! to prevent conflicts with the global memory.

  type t_storageBlock

    private

    ! Universally unique identifier
    type(t_uuid) :: ruuid

    ! An array of t_storageNode objects corresponding to the handles.
    ! Can be dynamically extended if there are not enough handles available.
    type(t_storageNode), dimension(:), pointer :: p_Rdescriptors => null()

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

    ! Maximum number of handles that were in use over the whole lifetime
    ! of this structure.
    integer :: nhandlesInUseMax = 0

    ! Maximum amount of memory that was in use over the whole lifetime
    ! of this structure.
    real(DP) :: dtotalMemMax = 0.0_DP

  end type t_storageBlock

  !</typeblock>

!</types>

!<globals>

  ! Global memory management structure
  type(t_storageBlock), private, save, target :: rbase

!</globals>

  interface storage_new
    module procedure storage_new1D
    module procedure storage_new1DFixed
    module procedure storage_new2D
    module procedure storage_new2DFixed
  end interface

  interface storage_realloc
    module procedure storage_realloc
    module procedure storage_reallocFixed
  end interface
  
  interface storage_getbase_int
    module procedure storage_getbase_int
    module procedure storage_getbase_intUBnd
    module procedure storage_getbase_intLUBnd
  end interface

  interface storage_getbase_single
    module procedure storage_getbase_single
    module procedure storage_getbase_singleUBnd
    module procedure storage_getbase_singleLUBnd
  end interface

  interface storage_getbase_double
    module procedure storage_getbase_double
    module procedure storage_getbase_doubleUBnd
    module procedure storage_getbase_doubleLUBnd
  end interface

  interface storage_getbase_logical
    module procedure storage_getbase_logical
    module procedure storage_getbase_logicalUBnd
    module procedure storage_getbase_logicalLUBnd
  end interface

  interface storage_getbase_char
    module procedure storage_getbase_char
    module procedure storage_getbase_charUBnd
    module procedure storage_getbase_charLUBnd
  end interface

  interface storage_getbase_int2D
    module procedure storage_getbase_int2D
    module procedure storage_getbase_int2DUBnd
    module procedure storage_getbase_int2DLUBnd
  end interface

  interface storage_getbase_single2D
    module procedure storage_getbase_single2D
    module procedure storage_getbase_single2DUBnd
    module procedure storage_getbase_single2DLUBnd
  end interface

  interface storage_getbase_double2D
    module procedure storage_getbase_double2D
    module procedure storage_getbase_double2DUBnd
    module procedure storage_getbase_double2DLUBnd
  end interface

  interface storage_getbase_logical2D
    module procedure storage_getbase_logical2D
    module procedure storage_getbase_logical2DUBnd
    module procedure storage_getbase_logical2DLUBnd
  end interface

  interface storage_getbase_char2D
    module procedure storage_getbase_char2D
    module procedure storage_getbase_char2DUBnd
    module procedure storage_getbase_char2DLUBnd
  end interface

  interface storage_getsize
    module procedure storage_getsize1D
    module procedure storage_getsize2D
  end interface

  interface storage_copy
    module procedure storage_copy
    module procedure storage_copy_explicit
    module procedure storage_copy_explicit2D
  end interface

contains

!************************************************************************

!<subroutine>

  subroutine storage_init(ihandleCount, ihandlesDelta, rheap)

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
  integer, intent(IN) :: ihandleCount

  ! OPTIONAL: Number of handles to increase the memory block by, if there are
  ! not enough handles available. Standard setting is 1/2*ihandleCount.
  integer, intent(IN), optional :: ihandlesDelta

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is initialised.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!</subroutine>

    ! local variables

    ! the real 'handle-delta'
    integer :: ihandles, ihDelta
  
    ! Pointer to the heap to initialise
    type(t_storageBlock), pointer :: p_rheap
    
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
      p_rheap => rbase
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
    
  end subroutine storage_init

!************************************************************************

!<subroutine>

  subroutine storage_done(rheap)

!<description>
  ! This routine cleans up the storage management. All data on the
  ! heap is released from memory.
!</description>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is cleaned up.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap to initialise
    type(t_storageBlock), pointer :: p_rheap
  
    integer :: i,ihandle
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    ! Delete all data from the heap
    do i = 1,size(p_rheap%p_Rdescriptors)
      ! Don't pass i as handle as storage_free will set the handle
      ! passed to it to 0!
      ihandle = i
      if (p_rheap%p_Rdescriptors(i)%idataType .ne. ST_NOHANDLE) &
          call storage_free(ihandle,rheap)
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
    
  end subroutine storage_done

!************************************************************************

!<function>

  integer function storage_newhandle (rheap) result(ihandle)

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
  type(t_storageBlock), intent(INOUT) :: rheap

!</inputoutput>

!</function>

    ! local variables
    type(t_storageNode), dimension(:), pointer :: p_Rdescriptors => null()
    integer, dimension(:), pointer :: p_IfreeHandles => null()
    integer :: i
    
    if (rheap%nhandlesTotal .le. 0) then
      call output_line ('Heap not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_newhandle')
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
    
  end function storage_newhandle

!************************************************************************

!<subroutine>

  subroutine storage_releasehandle (ihandle,rheap)

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
  type(t_storageBlock), intent(INOUT) :: rheap
!</inputoutput>

!</subroutine>

    type(t_storageNode), pointer :: p_rnode

    ! Where is the descriptor of the handle?
    p_rnode => rheap%p_Rdescriptors(ihandle)
    
    ! Subtract the memory amount from the statistics
    rheap%dtotalMem = rheap%dtotalMem - p_rnode%dmemBytes
    
    ! Clear the descriptor structure
    p_rnode%idataType = ST_NOHANDLE
    p_rnode%idimension = 0
    p_rnode%dmemBytes = 0.0_DP
    nullify(p_rnode%p_Fsingle1D)
    nullify(p_rnode%p_Ddouble1D)
    nullify(p_rnode%p_Iinteger1D)
    nullify(p_rnode%p_Blogical1D)
    nullify(p_rnode%p_Schar1D)
    nullify(p_rnode%p_Fsingle2D)
    nullify(p_rnode%p_Ddouble2D)
    nullify(p_rnode%p_Iinteger2D)
    nullify(p_rnode%p_Blogical2D)
    nullify(p_rnode%p_Schar2D)
    
    ! Handle ihandle is available now - put it to the list of available handles.
    rheap%p_ilastFreeHandle = mod(rheap%p_ilastFreeHandle,rheap%nhandlesTotal) + 1
    rheap%p_IfreeHandles (rheap%p_ilastFreeHandle) = ihandle
    
    rheap%ihandlesInUse = rheap%ihandlesInUse - 1

  end subroutine storage_releasehandle

!************************************************************************

!<subroutine>

  subroutine storage_initialiseNode (rstorageNode,cinitNewBlock,istartIndex,istopIndex)

!<description>
  ! Internal subroutine: Initialise the memory identified by storage
  ! node rstorageNode according to the constant cinitNewBlock.
!</description>

!<input>
  ! The storage node whose associated storage should be initialised.
  type(t_storageNode), intent(INOUT) :: rstorageNode

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO, ST_NEWBLOCK_NOINIT,
  ! ST_NEWBLOCK_ORDERED). Specifies how to initialise the data block associated
  ! to ihandle.
  integer, intent(IN) :: cinitNewBlock

  ! Start index from where to initialise; should usually be =1.
  ! For multidimensional arrays, this specifies the start index of the
  ! last dimension.
  integer(I32), intent(IN) :: istartIndex

  ! OPTIONAL: Stop index up to which to initialise; should usually be =SIZE(*).
  ! For multidimensional arrays, this specifies the stop index of the
  ! last dimension.
  integer(I32), intent(IN), optional :: istopIndex
!</input>

!</subroutine>

    ! variable for ordering 1,2,3,...,N
    integer :: iorder

    select case (rstorageNode%idimension)
    case (1)

      select case (cinitNewBlock)
      case (ST_NEWBLOCK_ZERO)
        ! Clear the vector if necessary
        if (present(istopIndex)) then
          select case (rstorageNode%idataType)
          case (ST_SINGLE)
            rstorageNode%p_Fsingle1D(istartIndex:istopIndex) = 0.0_SP
          case (ST_DOUBLE)
            rstorageNode%p_Ddouble1D(istartIndex:istopIndex) = 0.0_DP
          case (ST_INT)
            rstorageNode%p_Iinteger1D(istartIndex:istopIndex) = 0_I32
          case (ST_LOGICAL)
            rstorageNode%p_Blogical1D(istartIndex:istopIndex) = .false.
          case (ST_CHAR)
            rstorageNode%p_Schar1D(istartIndex:istopIndex) = achar(0)
          end select
        else
          select case (rstorageNode%idataType)
          case (ST_SINGLE)
            rstorageNode%p_Fsingle1D(istartIndex:) = 0.0_SP
          case (ST_DOUBLE)
            rstorageNode%p_Ddouble1D(istartIndex:) = 0.0_DP
          case (ST_INT)
            rstorageNode%p_Iinteger1D(istartIndex:) = 0_I32
          case (ST_LOGICAL)
            rstorageNode%p_Blogical1D(istartIndex:) = .false.
          case (ST_CHAR)
            rstorageNode%p_Schar1D(istartIndex:) = achar(0)
          end select
        end if

      case (ST_NEWBLOCK_ORDERED)
        ! Impose ordering 1,2,3,...,N if necessary
        if (present(istopIndex)) then
          select case (rstorageNode%idataType)
          case (ST_SINGLE)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Fsingle1D(iorder) = real(iorder,SP)
            end do
          case (ST_DOUBLE)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Ddouble1D(iorder) = real(iorder,DP)
            end do
          case (ST_INT)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Iinteger1D(iorder) = int(iorder,I32)
            end do
          case (ST_LOGICAL)
            call output_line ('Logical array can not be initialised with ST_NEWBLOCK_ORDERED!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'storage_initialiseNode')
            call sys_halt()
          case (ST_CHAR)
            call output_line ('Character array can not be initialised with ST_NEWBLOCK_ORDERED!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'storage_initialiseNode')
            call sys_halt()
          end select
        else
          select case (rstorageNode%idataType)
          case (ST_SINGLE)
            do iorder=istartIndex,ubound(rstorageNode%p_Fsingle1D,1)
              rstorageNode%p_Fsingle1D(iorder) = real(iorder,SP)
            end do
          case (ST_DOUBLE)
            do iorder=istartIndex,ubound(rstorageNode%p_Ddouble1D,1)
              rstorageNode%p_Ddouble1D(iorder) = real(iorder,DP)
            end do
          case (ST_INT)
            do iorder=istartIndex,ubound(rstorageNode%p_Iinteger1D,1)
              rstorageNode%p_Iinteger1D(iorder) = int(iorder,I32)
            end do
          case (ST_LOGICAL)
            call output_line ('Logical array can not be initialised with ST_NEWBLOCK_ORDERED!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'storage_initialiseNode')
            call sys_halt()
          case (ST_CHAR)
            call output_line ('Character array can not be initialised with ST_NEWBLOCK_ORDERED!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'storage_initialiseNode')
            call sys_halt()
          end select
        end if

      end select

    case (2)

      select case (cinitNewBlock)
      case (ST_NEWBLOCK_ZERO)
        ! Clear the vector
        if (present(istopIndex)) then
          select case (rstorageNode%idataType)
          case (ST_SINGLE)
            rstorageNode%p_Fsingle2D(:,istartIndex:istopIndex) = 0.0_SP
          case (ST_DOUBLE)
            rstorageNode%p_Ddouble2D(:,istartIndex:istopIndex) = 0.0_DP
          case (ST_INT)
            rstorageNode%p_Iinteger2D(:,istartIndex:istopIndex) = 0_I32
          case (ST_LOGICAL)
            rstorageNode%p_Blogical2D(:,istartIndex:istopIndex) = .false.
          case (ST_CHAR)
            rstorageNode%p_Schar2D(:,istartIndex:istopIndex) = achar(0)
          end select
        else
          select case (rstorageNode%idataType)
          case (ST_SINGLE)
            rstorageNode%p_Fsingle2D(:,istartIndex:) = 0.0_SP
          case (ST_DOUBLE)
            rstorageNode%p_Ddouble2D(:,istartIndex:) = 0.0_DP
          case (ST_INT)
            rstorageNode%p_Iinteger2D(:,istartIndex:) = 0_I32
          case (ST_LOGICAL)
            rstorageNode%p_Blogical2D(:,istartIndex:) = .false.
          case (ST_CHAR)
            rstorageNode%p_Schar2D(:,istartIndex:) = achar(0)
          end select
        end if

      case (ST_NEWBLOCK_ORDERED)
        call output_line ('Ordering not available for multidimensional array!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_initialiseNode')
        call sys_halt()

      end select

    case DEFAULT
      call output_line ('Unsupported dimension!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_initialiseNode')
      call sys_halt()

    end select

  end subroutine storage_initialiseNode

!************************************************************************

!<subroutine>

  subroutine storage_initialiseBlock (ihandle, cinitNewBlock,&
      rheap, istartIndex)

!<description>
  ! This routine initialises the memory associated to ihandle according
  ! to the constant cinitNewBlock.
!</description>

!<input>
  ! Handle of the memory block to initialise
  integer, intent(IN) :: ihandle

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO, ST_NEWBLOCK_NOINIT,
  ! ST_NEWBLOCK_ORDERED). Specifies how to initialise the data block associated
  ! to ihandle.
  integer, intent(IN) :: cinitNewBlock

  ! OPTIONAL: Start index of Block
  integer(I32), intent(IN), optional :: istartIndex
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap
!</inputoutput>

!</subroutine>

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    
    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Initialise
    if (present(istartIndex)) then
      call storage_initialiseNode (p_rnode,cinitNewBlock,istartIndex)
    else
      call storage_initialiseNode (p_rnode,cinitNewBlock,1_I32)
    end if

  end subroutine storage_initialiseBlock

!************************************************************************

!<subroutine>

  subroutine storage_new1D (scall, sname, isize, ctype, ihandle, &
                            cinitNewBlock, rheap)

!<description>
  !This routine reserves a 1D memory block of desired size and type.
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

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(OUT) :: ihandle

!</output>

!</subroutine>

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    character(LEN=SYS_NAMELEN) :: snameBackup
    
    if (isize .eq. 0) then
      call output_line ('isize=0!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'storage_new1D')
      ihandle = ST_NOHANDLE
      return
    end if

    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
  
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
    
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle1D(isize))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_SINGLE2BYTES)
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble1D(isize))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_DOUBLE2BYTES)
    case (ST_INT)
      allocate(p_rnode%p_Iinteger1D(isize))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_INT2BYTES)
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical1D(isize))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_LOGICAL2BYTES)
    case (ST_CHAR)
      allocate(p_rnode%p_Schar1D(isize))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_CHAR2BYTES)
    case DEFAULT
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new1D')
      call sys_halt()
    end select
    
    p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
    if (p_rheap%dtotalMem .gt. p_rheap%dtotalMemMax) &
        p_rheap%dtotalMemMax = p_rheap%dtotalMem
    
    ! Initialise the memory block
    call storage_initialiseBlock (ihandle, cinitNewBlock, rheap)
    
  end subroutine storage_new1D

!************************************************************************

!<subroutine>

  subroutine storage_new1Dfixed (scall, sname, ilbound, iubound,&
                            ctype, ihandle, cinitNewBlock, rheap)

!<description>
  !This routine reserves a 1D memory block of desired bounds and type.
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(IN) :: scall

  !clear name of data field
  character(LEN=*), intent(IN) :: sname

  !requested lower bound
  integer(I32), intent(IN) :: ilbound

  !requested upper bound
  integer(I32), intent(IN) :: iubound

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT,ST_NEWBLOCK_ORDERED)
  integer, intent(IN) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(OUT) :: ihandle

!</output>

!</subroutine>

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    integer(I32) :: isize
    character(LEN=SYS_NAMELEN) :: snameBackup
    
    isize=iubound-ilbound+1
    if (isize .eq. 0) then
      call output_line ('isize=0!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'storage_new1Dfixed')
      ihandle = ST_NOHANDLE
      return
    end if
    
    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
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
    
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle1D(ilbound:iubound))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_SINGLE2BYTES)
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble1D(ilbound:iubound))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_DOUBLE2BYTES)
    case (ST_INT)
      allocate(p_rnode%p_Iinteger1D(ilbound:iubound))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_INT2BYTES)
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical1D(ilbound:iubound))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_LOGICAL2BYTES)
    case (ST_CHAR)
      allocate(p_rnode%p_Schar1D(ilbound:iubound))
      p_rnode%dmemBytes = real(isize,DP)*real(ST_CHAR2BYTES)
    case DEFAULT
      call output_line ('Unsupported memory type!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'storage_new1Dfixed')
      call sys_halt()
    end select
    
    p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
    if (p_rheap%dtotalMem .gt. p_rheap%dtotalMemMax) &
        p_rheap%dtotalMemMax = p_rheap%dtotalMem
    
    ! Initialise the memory block
    call storage_initialiseBlock (ihandle, cinitNewBlock, rheap, ilbound)
    
  end subroutine storage_new1Dfixed

!************************************************************************

!<subroutine>

  subroutine storage_new2D (scall, sname, Isize, ctype, ihandle, &
                            cinitNewBlock, rheap)

!<description>
  !This routine reserves a 2D memory block of desired size and type.
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(IN) :: scall

  !clear name of data field
  character(LEN=*), intent(IN) :: sname

  !requested storage size for 1st and 2nd dimension
  integer(I32), dimension(2), intent(IN) :: Isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  integer, intent(IN) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(OUT) :: ihandle

!</output>

!</subroutine>

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    character(LEN=SYS_NAMELEN) :: snameBackup
    
    if ((Isize(1) .eq. 0) .or. (Isize(2) .eq. 0)) then
      call output_line ('Isize=0!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'storage_new2D')
      ihandle = 0
      return
    end if
    
    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
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
    
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle2D(Isize(1),Isize(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_SINGLE2BYTES)
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble2D(Isize(1),Isize(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_DOUBLE2BYTES)
    case (ST_INT)
      allocate(p_rnode%p_Iinteger2D(Isize(1),Isize(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_INT2BYTES)
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical2D(Isize(1),Isize(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_LOGICAL2BYTES)
    case (ST_CHAR)
      allocate(p_rnode%p_Schar2D(Isize(1),Isize(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_CHAR2BYTES)
    case DEFAULT
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new2D')
      call sys_halt()
    end select
    
    p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
    if (p_rheap%dtotalMem .gt. p_rheap%dtotalMemMax) &
        p_rheap%dtotalMemMax = p_rheap%dtotalMem
    
    ! Initialise the storage block
    call storage_initialiseBlock (ihandle, cinitNewBlock, rheap)
    
  end subroutine storage_new2D
  
!************************************************************************

!<subroutine>

  subroutine storage_new2Dfixed (scall, sname, Ilbound, Iubound, ctype,&
                            ihandle, cinitNewBlock, rheap)

!<description>
  !This routine reserves a 2D memory block of desired bounds and type.
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(IN) :: scall

  !clear name of data field
  character(LEN=*), intent(IN) :: sname

  !requested lower bounds for 1st and 2nd dimension
  integer(I32), dimension(2), intent(IN) :: Ilbound

  !requested upper bounds for 1st and 2nd dimension
  integer(I32), dimension(2), intent(IN) :: Iubound

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(IN) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  integer, intent(IN) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(OUT) :: ihandle

!</output>

!</subroutine>

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    integer(I32), dimension(2) :: Isize
    character(LEN=SYS_NAMELEN) :: snameBackup
    
    Isize=Iubound-Ilbound+1
    if ((Isize(1) .eq. 0) .or. (Isize(2) .eq. 0)) then
      call output_line ('Isize=0!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'storage_new2Dfixed')
      ihandle = 0
      return
    end if
    
    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
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
    
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_SINGLE2BYTES)
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_DOUBLE2BYTES)
    case (ST_INT)
      allocate(p_rnode%p_Iinteger2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_INT2BYTES)
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_LOGICAL2BYTES)
    case (ST_CHAR)
      allocate(p_rnode%p_Schar2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%dmemBytes = real(Isize(1),DP)*real(Isize(2),DP)*real(ST_CHAR2BYTES)
    case DEFAULT
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new2Dfixed')
      call sys_halt()
    end select
    
    p_rheap%dtotalMem = p_rheap%dtotalMem + p_rnode%dmemBytes
    if (p_rheap%dtotalMem .gt. p_rheap%dtotalMemMax) &
        p_rheap%dtotalMemMax = p_rheap%dtotalMem
    
    ! Initialise the storage block
    call storage_initialiseBlock (ihandle, cinitNewBlock, rheap, Ilbound(2))
    
  end subroutine storage_new2Dfixed

!************************************************************************

!<subroutine>

  subroutine storage_free (ihandle, rheap)

!<description>
  ! This routine releases a handle from a heap and deallocates the
  ! associated memory. ihandle is set to ST_NOHANDLE upon return.
!</description>

!<inputoutput>

  ! Handle of the memory block to be releases
  integer :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    
    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Releasing ST_NOHANDLE is not allowed!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_free')
      call sys_halt()
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Is the node associated at all?
    if (p_rnode%idataType .eq. ST_NOHANDLE) then
      call output_line ('Trying to release nonexistent handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_free')
      call output_line ('Handle number: '//trim(sys_siL(ihandle,11)), &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_free')
      call sys_halt()
    end if
    
    ! Release the memory assigned to that handle.
    if (associated(p_rnode%p_Fsingle1D))  deallocate(p_rnode%p_Fsingle1D)
    if (associated(p_rnode%p_Ddouble1D))  deallocate(p_rnode%p_Ddouble1D)
    if (associated(p_rnode%p_Iinteger1D)) deallocate(p_rnode%p_Iinteger1D)
    if (associated(p_rnode%p_Blogical1D)) deallocate(p_rnode%p_Blogical1D)
    if (associated(p_rnode%p_Schar1D))    deallocate(p_rnode%p_Schar1D)
    if (associated(p_rnode%p_Fsingle2D))  deallocate(p_rnode%p_Fsingle2D)
    if (associated(p_rnode%p_Ddouble2D))  deallocate(p_rnode%p_Ddouble2D)
    if (associated(p_rnode%p_Iinteger2D)) deallocate(p_rnode%p_Iinteger2D)
    if (associated(p_rnode%p_Blogical2D)) deallocate(p_rnode%p_Blogical2D)
    if (associated(p_rnode%p_Schar2D))    deallocate(p_rnode%p_Schar2D)
    
    ! Release the handle itself.
    call storage_releasehandle (ihandle,p_rheap)
    
    ! And finally reset the handle to ST_NOHANDLE.
    ihandle = ST_NOHANDLE
    
  end subroutine storage_free

!************************************************************************

!<subroutine>

  subroutine storage_clear (ihandle, rheap)

!<description>
  ! This routine clears an array identified by ihandle; all entries are
  ! overwritten by 0.
!</description>

!<inputoutput>

  ! Handle of the memory block to be releases
  integer :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    
    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'storage_clean')
      call sys_halt()
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Is the node associated at all?
    if (p_rnode%idataType .eq. ST_NOHANDLE) then
      call output_line ('Trying to clear nonexistent handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_clear')
      call output_line ('Handle number: '//trim(sys_siL(ihandle,11)), &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_clear')
      call sys_halt()
    end if
    
    ! What are we?
    select case (p_rnode%idimension)
    case (1)
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        p_rnode%p_Fsingle1D = 0.0_SP
      case (ST_DOUBLE)
        p_rnode%p_Ddouble1D = 0.0_DP
      case (ST_INT)
        p_rnode%p_Iinteger1D = 0_I32
      case (ST_LOGICAL)
        p_rnode%p_Blogical1D = .false.
      case (ST_CHAR)
        p_rnode%p_Schar1D = achar(0)
      end select
    case (2)
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        p_rnode%p_Fsingle2D = 0.0_SP
      case (ST_DOUBLE)
        p_rnode%p_Ddouble2D = 0.0_DP
      case (ST_INT)
        p_rnode%p_Iinteger2D = 0_I32
      case (ST_LOGICAL)
        p_rnode%p_Blogical2D = .false.
      case (ST_CHAR)
        p_rnode%p_Schar2D = achar(0)
      end select
    case DEFAULT
      call output_line ('Invalid dimension!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_clear')
      call sys_halt()
    end select
    
  end subroutine storage_clear

!************************************************************************

!<subroutine>

  subroutine storage_getsize1D (ihandle, isize, rheap)

!<description>
  ! Returns the length of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!<output>
  ! Length of the array identified by ihandle.
  integer(I32), intent(OUT) :: isize
!</output>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    
    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize1D')
      call sys_halt()
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Is the node associated at all?
    if (p_rnode%idataType .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize1D')
      call output_line ('Handle number: '//trim(sys_siL(ihandle,11)), &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize1D')
      call sys_halt()
    end if
    
    ! What are we?
    if (p_rnode%idimension .ne. 1) then
      call output_line ('Handle '//trim(sys_siL(ihandle,11))//' is not 1-dimensional!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize1D')
      call sys_halt()
    end if
    
    select case (p_rnode%idataType)
    case (ST_SINGLE)
      isize = size(p_rnode%p_Fsingle1D)
    case (ST_DOUBLE)
      isize = size(p_rnode%p_Ddouble1D)
    case (ST_INT)
      isize = size(p_rnode%p_Iinteger1D)
    case (ST_LOGICAL)
      isize = size(p_rnode%p_Blogical1D)
    case (ST_CHAR)
      isize = size(p_rnode%p_Schar1D)
    case DEFAULT
      call output_line ('Invalid data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize1D')
      call sys_halt()
    end select
    
  end subroutine storage_getsize1D

!************************************************************************

!<subroutine>

  subroutine storage_getsize2D (ihandle, isize, rheap)

!<description>
  ! Returns the length of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!<output>
  ! Length of each dimension of the array identified by ihandle.
  integer(I32), dimension(:), intent(OUT) :: Isize
!</output>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    
    ! Get the heap to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      call sys_halt()
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Is the node associated at all?
    if (p_rnode%idataType .eq. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      call output_line ('Handle number: '//trim(sys_siL(ihandle,11)), &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      call sys_halt()
    end if
    
    ! What are we?
    if (p_rnode%idimension .ne. 2) then
      call output_line ('Handle '//trim(sys_siL(ihandle,11))//' is not 2-dimensional!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize2D')
      call sys_halt()
    end if
    
    select case (p_rnode%idataType)
    case (ST_SINGLE)
      Isize = shape(p_rnode%p_Fsingle2D)
    case (ST_DOUBLE)
      Isize = shape(p_rnode%p_Ddouble2D)
    case (ST_INT)
      Isize = shape(p_rnode%p_Iinteger2D)
    case (ST_LOGICAL)
      Isize = shape(p_rnode%p_Blogical2D)
    case (ST_CHAR)
      Isize = shape(p_rnode%p_Schar2D)
    case DEFAULT
      call output_line ('Invalid data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize1D')
      call sys_halt()
    end select
    
  end subroutine storage_getsize2D

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int')
    call sys_halt()
  end if

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D

  end subroutine storage_getbase_int

!************************************************************************

!<subroutine>

  subroutine storage_getbase_intUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_intUBnd')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_intUBnd')
    call sys_halt()
  end if

  if ((ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D))) then
    call output_line ('Upper bound invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_intUBnd')
    call sys_halt()
  end if

  ! Get the pointer
  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D(:ubnd)

  end subroutine storage_getbase_intUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_intLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle
  
  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_intLUBnd')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_intLUBnd')
    call sys_halt()
  end if

  if ((lbnd .lt. 1) .or. &
      (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D)) .or. &
      (ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D)) .or. &
      (lbnd .gt. ubnd)) then
    call output_line ('Bound invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_intUBnd')
    call sys_halt()
  end if

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D(lbnd:ubnd)

  end subroutine storage_getbase_intLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single (ihandle, p_Sarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase__single')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single')
    call sys_halt()
  end if

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D

  end subroutine storage_getbase_single

!************************************************************************

!<subroutine>

  subroutine storage_getbase_singleUBnd (ihandle, p_Sarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_singleUBnd')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_singleUBnd')
    call sys_halt()
  end if

  if ((ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D))) then
    call output_line ('Upper bound invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_singleUBnd')
    call sys_halt()
  end if

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D(:ubnd)

  end subroutine storage_getbase_singleUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_singleLUBnd (ihandle, p_Sarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_singleLUBnd')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_singleLUBnd')
    call sys_halt()
  end if

  if ((lbnd .lt. 1) .or. &
      (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D)) .or. &
      (ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D)) .or. &
      (lbnd .gt. ubnd)) then
    call output_line ('Bound invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_singleUBnd')
    call sys_halt()
  end if

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D(lbnd:ubnd)

  end subroutine storage_getbase_singleLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double (ihandle, p_Darray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(DP), dimension(:), pointer :: p_Darray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double')
    call sys_halt()
  end if

  ! Get the pointer

  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D

  end subroutine storage_getbase_double

!************************************************************************

!<subroutine>

  subroutine storage_getbase_doubleUBnd (ihandle, p_Darray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(DP), dimension(:), pointer :: p_Darray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleUBnd')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleUBnd')
    call sys_halt()
  end if

  if ((ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D))) then
    call output_line ('Upper bound invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleUBnd')
    call sys_halt()
  end if

  ! Get the pointer
  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D(:ubnd)

  end subroutine storage_getbase_doubleUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_doubleLUBnd (ihandle, p_Darray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(DP), dimension(:), pointer :: p_Darray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleLUBnd')
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
    call output_line ('Wrong data format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleLUBnd')
    call sys_halt()
  end if

  if ((lbnd .lt. 1) .or. &
      (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D)) .or. &
      (ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D)) .or. &
      (lbnd .gt. ubnd)) then
    call output_line ('Bound invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleUBnd')
    call sys_halt()
  end if

  ! Get the pointer
  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D(lbnd:ubnd)

  end subroutine storage_getbase_doubleLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical (ihandle, p_Larray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:), pointer :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_logical: Wrong handle'
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
    print *,'storage_getbase_logical: Wrong data format!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D

  end subroutine storage_getbase_logical

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logicalUBnd (ihandle, p_Larray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:), pointer :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_logical: Wrong handle'
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
    print *,'storage_getbase_logical: Wrong data format!'
    call sys_halt()
  end if

  if ((ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D))) then
    print *,'storage_getbase_logicalUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D(:ubnd)

  end subroutine storage_getbase_logicalUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logicalLUBnd (ihandle, p_Larray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopts the given lower and upper bounds.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:), pointer :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_logical: Wrong handle'
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
    print *,'storage_getbase_logical: Wrong data format!'
    call sys_halt()
  end if

  if ((lbnd .lt. 1) .or. &
      (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .or. &
      (ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .or. &
      (lbnd .gt. ubnd)) then
    print *,'storage_getbase_logicalUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D(lbnd:ubnd)

  end subroutine storage_getbase_logicalLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char (ihandle, p_Carray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:), pointer :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_char: Wrong handle'
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
    print *,'storage_getbase_char: Wrong data format!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D

  end subroutine storage_getbase_char

!************************************************************************

!<subroutine>

  subroutine storage_getbase_charUBnd (ihandle, p_Carray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle
  
  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:), pointer :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_char: Wrong handle'
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
    print *,'storage_getbase_char: Wrong data format!'
    call sys_halt()
  end if

  if ((ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D))) then
    print *,'storage_getbase_charLUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D(:ubnd)

  end subroutine storage_getbase_charUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_charLUBnd (ihandle, p_Carray, lbnd, ubnd,rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopts the given lower and upper bounds.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:), pointer :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_char: Wrong handle'
    call sys_halt()
  end if

  if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
    print *,'storage_getbase_char: Wrong data format!'
    call sys_halt()
  end if

  if ((lbnd .lt. 1) .or. &
      (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .or. &
      (ubnd .lt. 1) .or. &
      (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .or. &
      (lbnd .gt. ubnd)) then
    print *,'storage_getbase_charLUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D(lbnd:ubnd)

  end subroutine storage_getbase_charLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int2D (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_int2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D

  end subroutine storage_getbase_int2D

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int2DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array and adopts the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_int2D: Wrong handle'
    call sys_halt()
  end if

  if ((ubnd .lt. 1) .or. &
      (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2))) then
    print *,'storage_getbase_int2DUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D(:,:ubnd)

  end subroutine storage_getbase_int2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int2DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! interger array and adopts the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_int2D: Wrong handle'
    call sys_halt()
  end if

  if ((lbnd .lt. 1) .or. &
      (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2)) .or. &
      (ubnd .lt. 1) .or. &
      (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2)) .or. &
      (lbnd .gt. ubnd)) then
    print *,'storage_getbase_int2DLUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D(:,lbnd:ubnd)

  end subroutine storage_getbase_int2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single2D (ihandle, p_Sarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! single precision array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_single2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D

  end subroutine storage_getbase_single2D

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single2DUBnd (ihandle, p_Sarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single precision array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_single2D: Wrong handle'
    call sys_halt()
  end if

  if ((ubnd .lt. 1) .or. &
      (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2))) then
    print *,'storage_getbase_single2DUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D(:,:ubnd)

  end subroutine storage_getbase_single2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single2DLUBnd (ihandle, p_Sarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! single precision array and adopt the given lower and upper bounds
  ! for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd
  
  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_single2D: Wrong handle'
    call sys_halt()
  end if

  if ((lbnd .lt. 1) .or. &
      (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2)) .or. &
      (ubnd .lt. 1) .or. &
      (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2)) .or. &
      (lbnd .gt. ubnd)) then
    print *,'storage_getbase_single2DLUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D(:,lbnd:ubnd)

  end subroutine storage_getbase_single2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double2D (ihandle, p_Darray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array.

!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_double2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D

  end subroutine storage_getbase_double2D

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double2DUBnd (ihandle, p_Darray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array and adopt the given upper bound for the second dimension.

!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_double2D: Wrong handle'
    call sys_halt()
  end if

  if ((ubnd .lt. 1) .or. &
      (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2))) then
    print *,'storage_getbase_double2DUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D(:,:ubnd)

  end subroutine storage_getbase_double2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double2DLUBnd (ihandle, p_Darray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array and adopt the given lower and upper bounds
  ! for the second dimension.

!</description>

!<input>
  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_double2D: Wrong handle'
    call sys_halt()
  end if

  if ((lbnd .lt. 1) .or. &
      (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2)) .or. &
      (ubnd .lt. 1) .or. &
      (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2)) .or. &
      (lbnd .gt. ubnd)) then
    print *,'storage_getbase_double2DLUBnd: Bounds invalid invalid!'
    call sys_halt()
  end if

  ! Get the pointer

  p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D(:,lbnd:ubnd)

  end subroutine storage_getbase_double2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical2D (ihandle, p_Larray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_logical2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D

  end subroutine storage_getbase_logical2D

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical2DUBnd (ihandle, p_Larray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_logical2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D(:,:ubnd)

  end subroutine storage_getbase_logical2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical2DLUBnd (ihandle, p_Larray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopt the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_logical2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D(:,lbnd:ubnd)

  end subroutine storage_getbase_logical2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char2D (ihandle, p_Carray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_char2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D

  end subroutine storage_getbase_char2D

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char2DUBnd (ihandle, p_Carray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_char2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D(:,:ubnd)

  end subroutine storage_getbase_char2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char2DLUBnd (ihandle, p_Carray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopt the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(IN) :: ihandle

  ! The lower bound
  integer, intent(IN) :: lbnd

  ! The upper bound
  integer, intent(IN) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .eq. ST_NOHANDLE) then
    print *,'storage_getbase_char2D: Wrong handle'
    call sys_halt()
  end if

  ! Get the pointer

  p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D(:,lbnd:ubnd)

  end subroutine storage_getbase_char2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_copy(h_source, h_dest, rheap)

!<description>
  ! This routine copies the information of one array to another.
  ! The structure of the arrays behind h_source and h_dest must be the
  ! same!
!</description>

!<input>
  ! Handle of the source array to copy
  integer, intent(IN) :: h_source

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  integer, intent(INOUT) :: h_dest
!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rsource, p_rdest
  integer(I32) :: i,j
  integer(I32), dimension(2) :: Isize

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (h_source .eq. ST_NOHANDLE) then
      print *,'storage_copy: Wrong handle'
      call sys_halt()
    end if
    if (.not. associated(p_rheap%p_Rdescriptors)) then
      print *,'storage_copy: Heap not initialised!'
      call sys_halt()
    end if

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    if (h_dest .eq. ST_NOHANDLE) then
      ! Create a new array in the same size and structure
      ! as h_source.
      select case (p_rsource%idimension)
      case (1)
        select case (p_rsource%idataType)
        case (ST_DOUBLE)
          call storage_new ('storage_copy',p_rsource%sname,&
                            int(size(p_rsource%p_Ddouble1D),I32),&
                            ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_SINGLE)
          call storage_new ('storage_copy',p_rsource%sname,&
                            int(size(p_rsource%p_Fsingle1D),I32),&
                            ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT)
          call storage_new ('storage_copy',p_rsource%sname,&
                            int(size(p_rsource%p_Iinteger1D),I32),&
                            ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_LOGICAL)
          call storage_new ('storage_copy',p_rsource%sname,&
                            int(size(p_rsource%p_Blogical1D),I32),&
                            ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_CHAR)
          call storage_new ('storage_copy',p_rsource%sname,&
                            int(size(p_rsource%p_Schar1D),I32),&
                            ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        end select
      case (2)
        select case (p_rsource%IdataType)
        case (ST_DOUBLE)
          Isize = ubound(p_rsource%p_Ddouble2D)   ! =SIZE(...) here
          call storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_SINGLE)
          Isize = ubound(p_rsource%p_Fsingle2D)   ! =SIZE(...) here
          call storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT)
          Isize = ubound(p_rsource%p_Iinteger2D)      ! =SIZE(...) here
          call storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_LOGICAL)
          Isize = ubound(p_rsource%p_Blogical2D)      ! =SIZE(...) here
          call storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_CHAR)
          Isize = ubound(p_rsource%p_Schar2D)      ! =SIZE(...) here
          call storage_new ('storage_copy', p_rsource%sname, Isize,&
                            ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        end select
      end select

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it's correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      
    end if

    p_rdest => p_rheap%p_Rdescriptors(h_dest)

    ! 1D/2D the same?
    if (p_rsource%idimension .ne. p_rdest%idimension) then
      print *,'storage_copy: Structure different!'
      call sys_halt()
    end if

    ! What is to copy
    select case (p_rsource%idimension)
    case (1)
      select case (p_rsource%idataType)
      case (ST_DOUBLE)
        select case (p_rdest%idataType)
        case (ST_DOUBLE)
          call lalg_copyVectorDble (p_rsource%p_Ddouble1D,p_rdest%p_Ddouble1D)
        case (ST_SINGLE)
          call lalg_copyVectorDblSngl (p_rsource%p_Ddouble1D,p_rdest%p_Fsingle1D)
        case DEFAULT
          print *,'storage_copy: Unsupported data type combination'
          call sys_halt()
        end select

      case (ST_SINGLE)
        select case (p_rdest%idataType)
        case (ST_DOUBLE)
          call lalg_copyVectorSnglDbl (p_rsource%p_Fsingle1D,p_rdest%p_Ddouble1D)
        case (ST_SINGLE)
          call lalg_copyVectorSngl (p_rsource%p_Fsingle1D,p_rdest%p_Fsingle1D)
        case DEFAULT
          print *,'storage_copy: Unsupported data type combination'
          call sys_halt()
        end select

      case (ST_INT)
        if (p_rdest%idataType .eq. ST_INT) then
          call lalg_copyVectorInt (p_rsource%p_Iinteger1D,p_rdest%p_Iinteger1D)
        else
          print *,'storage_copy: Unsupported data type combination'
          call sys_halt()
        end if

      case (ST_LOGICAL)
        if (p_rdest%idataType .eq. ST_LOGICAL) then
          ! Copy by hand
          do i=1,size(p_rsource%p_Blogical1D)
            p_rdest%p_Blogical1D(i) = p_rsource%p_Blogical1D(i)
          end do
        else
          print *,'storage_copy: Unsupported data type combination'
          call sys_halt()
        end if

      case (ST_CHAR)
        if (p_rdest%idataType .eq. ST_CHAR) then
          ! Copy by hand
          do i=1,size(p_rsource%p_Schar1D)
            p_rdest%p_Schar1D(i) = p_rsource%p_Schar1D(i)
          end do
        else
          print *,'storage_copy: Unsupported data type combination'
          call sys_halt()
        end if

      case DEFAULT
        print *,'storage_copy: Unknown data type'
        call sys_halt()
      end select

    case (2)
      select case (p_rsource%idataType)
      case (ST_DOUBLE)
        if ((ubound(p_rsource%p_Ddouble2D,1) .ne. ubound(p_rdest%p_Ddouble2D,1)) .or.&
            (ubound(p_rsource%p_Ddouble2D,2) .ne. ubound(p_rdest%p_Ddouble2D,2))) then
          print *,'storage_copy: Structure different!'
          call sys_halt()
        end if

        select case (p_rdest%idataType)
        case (ST_DOUBLE)
          ! Copy by hand
          do j=1,ubound(p_rsource%p_Ddouble2D,2)
            do i=1,ubound(p_rsource%p_Ddouble2D,1)
              p_rdest%p_Ddouble2D(i,j) = p_rsource%p_Ddouble2D(i,j)
            end do
          end do

        case (ST_SINGLE)
          ! Copy by hand
          do j=1,ubound(p_rsource%p_Fsingle2D,2)
            do i=1,ubound(p_rsource%p_Fsingle2D,1)
              p_rdest%p_Fsingle2D(i,j) = p_rsource%p_Ddouble2D(i,j)
            end do
          end do

        case (ST_INT)
          ! Copy by hand
          do j=1,ubound(p_rsource%p_Iinteger2D,2)
            do i=1,ubound(p_rsource%p_Iinteger2D,1)
              p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Ddouble2D(i,j)
            end do
          end do

        case DEFAULT ! Might be ST_LOGICAL or ST_CHAR
          print *,'storage_copy: Unsupported data type combination'
          call sys_halt()

        end select

      case (ST_SINGLE)
        if ((ubound(p_rsource%p_Fsingle2D,1) .ne. ubound(p_rdest%p_Fsingle2D,1)) .or.&
            (ubound(p_rsource%p_Fsingle2D,2) .ne. ubound(p_rdest%p_Fsingle2D,2))) then
          print *,'storage_copy: Structure different!'
          call sys_halt()
        end if

        select case (p_rdest%idataType)
        case (ST_DOUBLE)
          ! Copy by hand
          do j=1,ubound(p_rsource%p_Ddouble2D,2)
            do i=1,ubound(p_rsource%p_Ddouble2D,1)
              p_rdest%p_Ddouble2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            end do
          end do

        case (ST_SINGLE)
          ! Copy by hand
          do j=1,ubound(p_rsource%p_Fsingle2D,2)
            do i=1,ubound(p_rsource%p_Fsingle2D,1)
              p_rdest%p_Fsingle2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            end do
          end do

        case (ST_INT)
          ! Copy by hand
          do j=1,ubound(p_rsource%p_Iinteger2D,2)
            do i=1,ubound(p_rsource%p_Iinteger2D,1)
              p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Fsingle2D(i,j)
            end do
          end do

        ! Might be ST_LOGICAL or ST_CHAR
        case DEFAULT
          print *,'storage_copy: Unsupported data type combination'
          call sys_halt()

        end select

      case (ST_INT)
        if ((ubound(p_rsource%p_Iinteger2D,1) .ne. ubound(p_rdest%p_Iinteger2D,1)) .or.&
            (ubound(p_rsource%p_Iinteger2D,2) .ne. ubound(p_rdest%p_Iinteger2D,2))) then
          print *,'storage_copy: Structure different!'
          call sys_halt()
        end if
        if (p_rdest%idataType .ne. ST_INT) then
          print *,'storage_copy: unsupported data type combination'
          call sys_halt()
        end if

        ! Copy by hand
        do j=1,ubound(p_rsource%p_Iinteger2D,2)
          do i=1,ubound(p_rsource%p_Iinteger2D,1)
            p_rdest%p_Iinteger2D(i,j) = p_rsource%p_Iinteger2D(i,j)
          end do
        end do

      case (ST_LOGICAL)
        if ((ubound(p_rsource%p_Blogical2D,1) .ne. ubound(p_rdest%p_Blogical2D,1)) .or.&
            (ubound(p_rsource%p_Blogical2D,2) .ne. ubound(p_rdest%p_Blogical2D,2))) then
          print *,'storage_copy: Structure different!'
          call sys_halt()
        end if
        if (p_rdest%idataType .ne. ST_LOGICAL) then
          print *,'storage_copy: unsupported data type combination'
          call sys_halt()
        end if

        ! Copy by hand
        do j=1,ubound(p_rsource%p_Blogical2D,2)
          do i=1,ubound(p_rsource%p_Blogical2D,1)
            p_rdest%p_Blogical2D(i,j) = p_rsource%p_Blogical2D(i,j)
          end do
        end do

      case (ST_CHAR)
        if ((ubound(p_rsource%p_Schar2D,1) .ne. ubound(p_rdest%p_Schar2D,1)) .or.&
            (ubound(p_rsource%p_Schar2D,2) .ne. ubound(p_rdest%p_Schar2D,2))) then
          print *,'storage_copy: Structure different!'
          call sys_halt()
        end if
        if (p_rdest%idataType .ne. ST_CHAR) then
          print *,'storage_copy: unsupported data type combination'
          call sys_halt()
        end if

        ! Copy by hand
        do j=1,ubound(p_rsource%p_Schar2D,2)
          do i=1,ubound(p_rsource%p_Schar2D,1)
            p_rdest%p_Schar2D(i,j) = p_rsource%p_Schar2D(i,j)
          end do
        end do

      case DEFAULT
        print *,'storage_copy: Unknown data type'
        call sys_halt()
      end select
    end select

  end subroutine storage_copy

!************************************************************************

!<subroutine>

  subroutine storage_copy_explicit(h_source, h_dest, istart_source, &
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
  integer, intent(IN) :: h_source

  ! First entry of the source array to copy
  integer, intent(IN) :: istart_source

  ! First entry of the destination array where to copy
  integer, intent(IN) :: istart_dest

  ! Length of the array to copy
  integer, intent(IN) :: ilength

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  integer, intent(INOUT) :: h_dest
!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rsource, p_rdest
    integer(I32) :: i
    character(LEN=SYS_NAMELEN) :: sname = ''

    ! Check if the start address is positive
    if (istart_source <= 0 .or. istart_dest <= 0) then
      print *, 'storage_copy_explicit: start address must be positive'
      call sys_halt()
    end if

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (h_source .eq. ST_NOHANDLE) then
      print *,'storage_copy_explicit: Wrong handle'
      call sys_halt()
    end if
    if (.not. associated(p_rheap%p_Rdescriptors)) then
      print *,'storage_copy_explicit: Heap not initialised!'
      call sys_halt()
    end if

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    if (h_dest .eq. ST_NOHANDLE) then
      ! Create a new array in the same size and structure
      ! as h_source.
      if (p_rsource%idimension /= 1) then
        print *, 'storage_copy_explicit: only 1D arrays are allowed'
        call sys_halt()
      end if

      select case (p_rsource%idataType)
      case (ST_DOUBLE)
         call storage_new ('storage_copy_explicit',p_rsource%sname,&
              int(size(p_rsource%p_Ddouble1D),I32),&
              ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_SINGLE)
         call storage_new ('storage_copy_explicit',p_rsource%sname,&
              int(size(p_rsource%p_Fsingle1D),I32),&
              ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT)
         call storage_new ('storage_copy_explicit',p_rsource%sname,&
              int(size(p_rsource%p_Iinteger1D),I32),&
              ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_LOGICAL)
         call storage_new ('storage_copy_explicit',p_rsource%sname,&
              int(size(p_rsource%p_Blogical1D),I32),&
              ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_CHAR)
         call storage_new ('storage_copy_explicit',p_rsource%sname,&
              int(size(p_rsource%p_Schar1D),I32),&
              ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      end select

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it's correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      
    end if

    p_rdest => p_rheap%p_Rdescriptors(h_dest)

    ! 1D/2D the same?
    if (p_rsource%idimension .ne. p_rdest%idimension) then
      print *,'storage_copy_explicit: Structure different!'
      call sys_halt()
    end if

    ! What is to copy
    select case (p_rsource%idataType)
    case (ST_DOUBLE)
       select case (p_rdest%idataType)
       case (ST_DOUBLE)
          if (istart_source+ilength-1 > size(p_rsource%p_Ddouble1D) .or. &
               istart_dest+ilength-1 > size(p_rdest%p_Ddouble1D)) then
            print *, 'storage_copy_explicit: Subarrays incompatible!'
            call sys_halt()
          end if
          ! Copy by hand
          do i=1,ilength
            p_rdest%p_Ddouble1D(istart_dest+i-1) = &
              p_rsource%p_Ddouble1D(istart_source+i-1)
          end do
       case (ST_SINGLE)
          if (istart_source+ilength-1 > size(p_rsource%p_Ddouble1D) .or. &
               istart_dest+ilength-1 > size(p_rdest%p_Fsingle1D)) then
            print *, 'storage_copy_explicit: Subarrays incompatible!'
            call sys_halt()
          end if
          ! Copy by hand
          do i=1,ilength
            p_rdest%p_Fsingle1D(istart_dest+i-1) = &
              p_rsource%p_Ddouble1D(istart_source+i-1)
          end do
       case DEFAULT
          print *,'storage_copy_explicit: Unsupported data type combination'
          call sys_halt()
       end select

    case (ST_SINGLE)
       select case (p_rdest%idataType)
       case (ST_DOUBLE)
          if (istart_source+ilength-1 > size(p_rsource%p_Fsingle1D) .or. &
               istart_dest+ilength-1 > size(p_rdest%p_Ddouble1D)) then
             print *, 'storage_copy_explicit: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do i=1,ilength
             p_rdest%p_Ddouble1D(istart_dest+i-1) = &
               p_rsource%p_Fsingle1D(istart_source+i-1)
          end do
       case (ST_SINGLE)
          if (istart_source+ilength-1 > size(p_rsource%p_Fsingle1D) .or. &
               istart_dest+ilength-1 > size(p_rdest%p_Fsingle1D)) then
             print *, 'storage_copy_explicit: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do i=1,ilength
             p_rdest%p_Fsingle1D(istart_dest+i-1) = &
               p_rsource%p_Fsingle1D(istart_source+i-1)
          end do
       case DEFAULT
          print *,'storage_copy_explicit: Unsupported data type combination'
          call sys_halt()
       end select

    case (ST_INT)
       if (p_rdest%idataType .eq. ST_INT) then
          if (istart_source+ilength-1 > size(p_rsource%p_Iinteger1D) .or. &
               istart_dest+ilength-1 > size(p_rdest%p_Iinteger1D)) then
             print *, 'storage_copy_explicit: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do i=1,ilength
             p_rdest%p_Iinteger1D(istart_dest+i-1) = p_rsource%p_Iinteger1D(istart_source+i-1)
          end do
       else
          print *,'storage_copy_explicit: Unsupported data type combination'
          call sys_halt()
       end if

    case (ST_LOGICAL)
       if (p_rdest%idataType .eq. ST_LOGICAL) then
          if (istart_source+ilength-1 > size(p_rsource%p_Blogical1D) .or. &
               istart_dest+ilength-1 > size(p_rdest%p_Blogical1D)) then
             print *, 'storage_copy_explicit: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do i=1,ilength
             p_rdest%p_Blogical1D(istart_dest+i-1) = p_rsource%p_Blogical1D(istart_source+i-1)
          end do
       else
          print *,'storage_copy_explicit: Unsupported data type combination'
          call sys_halt()
       end if

    case (ST_CHAR)
       if (p_rdest%idataType .eq. ST_CHAR) then
          if (istart_source+ilength-1 > size(p_rsource%p_Schar1D) .or. &
               istart_dest+ilength-1 > size(p_rdest%p_Schar1D)) then
             print *, 'storage_copy_explicit: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do i=1,ilength
             p_rdest%p_Schar1D(istart_dest+i-1) = p_rsource%p_Schar1D(istart_source+i-1)
          end do
       else
          print *,'storage_copy_explicit: Unsupported data type combination'
          call sys_halt()
       end if

    case DEFAULT
       print *,'storage_copy_explicit: Unknown data type'
       call sys_halt()
    end select

  end subroutine storage_copy_explicit

!************************************************************************

!<subroutine>

  subroutine storage_copy_explicit2D(h_source, h_dest, Istart_source, &
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
  integer, intent(IN) :: h_source

  ! First entry of the source array to copy
  integer, dimension(2), intent(IN) :: Istart_source

  ! First entry of the destination array where to copy
  integer, dimension(2), intent(IN) :: Istart_dest

  ! Length of the array to copy
  integer, dimension(2), intent(IN) :: Ilength

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  integer, intent(INOUT) :: h_dest
!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rsource, p_rdest
  integer(I32) :: i,j
  integer(I32), dimension(2) :: Isize

    ! Check if the start address is positive
    if (any(istart_source <= 0) .or. any(istart_dest <= 0)) then
      print *, 'storage_copy_explicit: start address must be positive'
      call sys_halt()
    end if

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (h_source .eq. ST_NOHANDLE) then
      print *,'storage_copy_explicit: Wrong handle'
      call sys_halt()
    end if
    if (.not. associated(p_rheap%p_Rdescriptors)) then
      print *,'storage_copy_explicit: Heap not initialised!'
      call sys_halt()
    end if

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    if (h_dest .eq. ST_NOHANDLE) then
      ! Create a new array in the same size and structure
      ! as h_source.
      if (p_rsource%idimension /= 2) then
        print *, 'storage_copy_explicit: only 1D arrays are allowed'
        call sys_halt()
      end if

      select case (p_rsource%IdataType)
      case (ST_DOUBLE)
         Isize = ubound(p_rsource%p_Ddouble2D)   ! =SIZE(...) here
         call storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_SINGLE)
         Isize = ubound(p_rsource%p_Fsingle2D)   ! =SIZE(...) here
         call storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT)
         Isize = ubound(p_rsource%p_Iinteger2D)      ! =SIZE(...) here
         call storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_LOGICAL)
         Isize = ubound(p_rsource%p_Blogical2D)      ! =SIZE(...) here
         call storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_CHAR)
         Isize = ubound(p_rsource%p_Schar2D)      ! =SIZE(...) here
         call storage_new ('storage_copy', p_rsource%sname, Isize,&
              ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      end select

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it's correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      
    end if

    p_rdest => p_rheap%p_Rdescriptors(h_dest)

    ! 1D/2D the same?
    if (p_rsource%idimension .ne. p_rdest%idimension) then
      print *,'storage_copy_explicit: Structure different!'
      call sys_halt()
    end if

    ! What is to copy
    select case (p_rsource%idataType)
    case (ST_DOUBLE)
       select case (p_rdest%idataType)
       case (ST_DOUBLE)
          if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
               Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
               Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Ddouble2D,1) .or. &
               Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Ddouble2D,2)) then
             print *, 'storage_copy_explicit2D: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do j=1,Ilength(2)
             do i=1,Ilength(1)
                p_rdest%p_Ddouble2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             end do
          end do

       case (ST_SINGLE)
          if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
               Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
               Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Fsingle2D,1) .or. &
               Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Fsingle2D,2)) then
             print *, 'storage_copy_explicit2D: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do j=1,Ilength(2)
             do i=1,Ilength(1)
                p_rdest%p_Fsingle2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             end do
          end do

       case (ST_INT)
          if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
               Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
               Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Iinteger2D,1) .or. &
               Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Iinteger2D,2)) then
             print *, 'storage_copy_explicit2D: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do j=1,Ilength(2)
             do i=1,Ilength(1)
                p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             end do
          end do

       case DEFAULT ! Might be ST_LOGICAL or ST_CHAR
          print *,'storage_copy_explicit2D: Unsupported data type combination'
          call sys_halt()

       end select

    case (ST_SINGLE)
       select case (p_rdest%idataType)
       case (ST_DOUBLE)
          if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
               Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
               Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Ddouble2D,1) .or. &
               Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Ddouble2D,2)) then
             print *, 'storage_copy_explicit2D: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do j=1,Ilength(2)
             do i=1,Ilength(1)
                p_rdest%p_Ddouble2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             end do
          end do

       case (ST_SINGLE)
          if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
               Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
               Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Fsingle2D,1) .or. &
               Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Fsingle2D,2)) then
             print *, 'storage_copy_explicit2D: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do j=1,Ilength(2)
             do i=1,Ilength(1)
                p_rdest%p_Fsingle2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             end do
          end do

       case (ST_INT)
          if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
               Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
               Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Iinteger2D,1) .or. &
               Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Iinteger2D,2)) then
             print *, 'storage_copy_explicit2D: Subarrays incompatible!'
             call sys_halt()
          end if
          ! Copy by hand
          do j=1,Ilength(2)
             do i=1,Ilength(1)
                p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                     p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
             end do
          end do

       case DEFAULT ! Might be ST_LOGICAL or ST_CHAR
          print *,'storage_copy_explicit2D: Unsupported data type combination'
          call sys_halt()

       end select

    case (ST_INT)
       if (p_rdest%idataType .ne. ST_INT) then
          print *,'storage_copy_explicit2D: unsupported data type combination'
          call sys_halt()
       end if
       if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Iinteger2D,2)) then
          print *, 'storage_copy_explicit2D: Subarrays incompatible!'
          call sys_halt()
       end if

       ! Copy by hand
       do j=1,Ilength(2)
          do i=1,Ilength(1)
             p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                  p_rsource%p_Iinteger2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
       end do

    case (ST_LOGICAL)
       if (p_rdest%idataType .ne. ST_LOGICAL) then
          print *,'storage_copy_explicit2D: unsupported data type combination'
          call sys_halt()
       end if
       if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Blogical2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Blogical2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Blogical2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Blogical2D,2)) then
          print *, 'storage_copy_explicit2D: Subarrays incompatible!'
          call sys_halt()
       end if

       ! Copy by hand
       do j=1,Ilength(2)
          do i=1,Ilength(1)
             p_rdest%p_Blogical2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                  p_rsource%p_Blogical2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
       end do

    case (ST_CHAR)
       if (p_rdest%idataType .ne. ST_CHAR) then
          print *,'storage_copy_explicit2D: unsupported data type combination'
          call sys_halt()
       end if
       if (Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Schar2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Schar2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1 > ubound(p_rdest%p_Schar2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1 > ubound(p_rdest%p_Schar2D,2)) then
          print *, 'storage_copy_explicit2D: Subarrays incompatible!'
          call sys_halt()
       end if

       ! Copy by hand
       do j=1,Ilength(2)
          do i=1,Ilength(1)
             p_rdest%p_Schar2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                  p_rsource%p_Schar2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
       end do

    case DEFAULT
       print *,'storage_copy: Unknown data type'
       call sys_halt()
    end select

  end subroutine storage_copy_explicit2D

!************************************************************************

!<subroutine>

  subroutine storage_info(bprintHandles,rheap)

!<description>
  ! This routine prints information about the current memory consumption
  ! in a memory block to screen.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the handles still remaining in the
  ! heap together with their names are printed to the terminal.
  logical, intent(IN), optional :: bprintHandles

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!</subroutine>

  ! local variables
  integer :: i

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    call output_line ('Heap statistics:')
    call output_line ('----------------')
    if (present(bprintHandles)) then
      if (bprintHandles .and. (p_rheap%ihandlesInUse .gt. 0)) then
        call output_line ('Handles on the heap: ')
        call output_lbrk ()
        ! Loop through the heap and search allocated handles
        do i=1,size(p_rheap%p_IfreeHandles)
          if (p_rheap%p_Rdescriptors(i)%idataType .ne. ST_NOHANDLE) then
            if (p_rheap%p_Rdescriptors(i)%idimension .eq. 1) then
              call output_line ( &
                   'Handle ' // trim(sys_siL(i,10)) // ', 1D, Length=' // &
                   trim(sys_siL(int(p_rheap%p_Rdescriptors(i)%dmemBytes),15)) //&
                   ', Type=' // trim(sys_siL(p_rheap%p_Rdescriptors(i)%idataType,15)) //&
                   ' Name=' // trim(adjustl(p_rheap%p_Rdescriptors(i)%sname)) )
            else
              call output_line ( &
                   'Handle ' // trim(sys_siL(i,10)) // ', 2D, Length=' // &
                   trim(sys_siL(int(p_rheap%p_Rdescriptors(i)%dmemBytes),15)) // &
                   ', Type=' // trim(sys_siL(p_rheap%p_Rdescriptors(i)%idataType,15)) //&
                   ' Name=' // trim(adjustl(p_rheap%p_Rdescriptors(i)%sname)) )
            end if
          end if
        end do
        call output_lbrk ()
      end if
    end if

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
  end subroutine storage_info

!************************************************************************

!<subroutine>

  subroutine storage_getdatatype (ihandle, idatatype, rheap)

!<description>
  ! Returns the datatype of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!<output>
  ! Datatype of the array identified by ihandle.
  integer(I32), intent(OUT) :: idatatype
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .le. ST_NOHANDLE) then
    call output_line ('Handle invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getdatatype')
    call sys_halt()
  end if

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  select case (p_rnode%idataType)
  case (ST_SINGLE)
    idatatype = ST_SINGLE
  case (ST_DOUBLE)
    idatatype = ST_DOUBLE
  case (ST_INT)
    idatatype = ST_INT
  case (ST_LOGICAL)
    idatatype = ST_LOGICAL
  case (ST_CHAR)
    idatatype = ST_CHAR
  case (ST_NOHANDLE)
    call output_line ('Handle invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getdatatype')
    call output_line ('Handle number: '//trim(sys_siL(ihandle,11)), &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getdatatype')
    call sys_halt()
  case DEFAULT
    call output_line ('Invalid data type!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getdatatype')
    call sys_halt()
  end select

  end subroutine storage_getdatatype

!************************************************************************

!<subroutine>

  subroutine storage_getdimension (ihandle, idimension, rheap)

!<description>
  ! Returns the dimension of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!<output>
  ! Dimension of the array identified by ihandle.
  integer(I32), intent(OUT) :: idimension
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .le. ST_NOHANDLE) then
    call output_line ('Handle invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getdimension')
    call sys_halt()
  end if

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)
  
  idimension = p_rnode%idimension

  end subroutine storage_getdimension

!************************************************************************

!<subroutine>

  subroutine storage_getname (ihandle, sname, rheap)

!<description>
  ! Returns the name of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(IN) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap
!</input>

!<output>
  ! Name of the array identified by ihandle.
  character(LEN=*), intent(OUT) :: sname
!</output>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode

  ! Get the heap to use - local or global one.

  if(present(rheap)) then
    p_rheap => rheap
  else
    p_rheap => rbase
  end if

  if (ihandle .le. ST_NOHANDLE) then
    call output_line ('Handle invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getname')
    call sys_halt()
  end if

  ! Where is the descriptor of the handle?
  p_rnode => p_rheap%p_Rdescriptors(ihandle)

  sname = p_rnode%sname

  end subroutine storage_getname

!************************************************************************

!<subroutine>

  subroutine storage_realloc (scall, isize, ihandle, cinitNewBlock, bcopy, rheap)

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
  character(LEN=*), intent(IN) :: scall

  ! Requested storage size for the memory block / the new size of the last
  ! dimension in the memory block identified by ihandle
  integer(I32), intent(IN) :: isize

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO,
  ! ST_NEWBLOCK_NOINIT, ST_NEWBLOCK_ORDERED).
  ! Specifies how to initialise memory if isize > original array size.
  integer, intent(IN) :: cinitNewBlock

  ! OPTIONAL: Copy old data.
  ! =TRUE: Copy data of old array to the new one.
  ! =FALSE: Reallocate memory, don't copy old data.
  ! If not specified, TRUE is assumed.
  logical, intent(IN), optional :: bcopy

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

  ! Handle of the memory block.
  integer, intent(INOUT) :: ihandle

!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode

    ! New storage node
    type(t_storageNode) :: rstorageNode

    ! size of the old 1-dimensional array
    integer(I32) :: isizeOld

    ! size of the 1-dimensional array to be copied
    integer :: isizeCopy

    ! size of the old 2-dimensional array
    integer(I32), dimension(2) :: Isize2Dold

    integer(I32) :: i,j

    logical :: bcopyData

    if (isize .eq. 0) then
      ! Ok, not much to do...
      call storage_free(ihandle,rheap)
      return
    end if

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    ! Copy old data?

    if (present(bcopy)) then
      bcopyData = bcopy
    else
      bcopyData = .true.
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Copy the data of the old storage node to rstorageNode. That way
    ! we prepare a new storage node and will replace the old.
    rstorageNode = p_rnode

    ! Are we 1D or 2D?
    select case(p_rnode%idimension)

    case (1)

      ! Get the size of the old storage node.
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        isizeOld = size(p_rnode%p_Fsingle1D)
      case (ST_DOUBLE)
        isizeOld = size(p_rnode%p_Ddouble1D)
      case (ST_INT)
        isizeOld = size(p_rnode%p_Iinteger1D)
      case (ST_LOGICAL)
        isizeOld = size(p_rnode%p_Blogical1D)
      case (ST_CHAR)
        isizeOld = size(p_rnode%p_Schar1D)
      end select

      ! Do we really have to change anything?
      if (isize == isizeOld) return

      ! Allocate new memory and initialise it - if it's larger than the old
      ! memory block.

      select case (rstorageNode%idataType)
      case (ST_SINGLE)
        allocate(rstorageNode%p_Fsingle1D(isize))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_SINGLE2BYTES,DP)
      case (ST_DOUBLE)
        allocate(rstorageNode%p_Ddouble1D(isize))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_DOUBLE2BYTES,DP)
      case (ST_INT)
        allocate(rstorageNode%p_Iinteger1D(isize))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_INT2BYTES,DP)
      case (ST_LOGICAL)
        allocate(rstorageNode%p_Blogical1D(isize))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_LOGICAL2BYTES,DP)
      case (ST_CHAR)
        allocate(rstorageNode%p_Schar1D(isize))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_CHAR2BYTES,DP)
      case DEFAULT
        print *,'Error: unknown mem type'
        call sys_halt()
      end select

      if (isize > isizeOld) &
        call storage_initialiseNode (rstorageNode,cinitNewBlock,isizeOld+1_I32)

      ! Copy old data?
      if (bcopyData) then
        isizeCopy=min(isize,isizeOld)
        select case (rstorageNode%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorSngl (p_rnode%p_Fsingle1D(1:isizeCopy),&
                                    rstorageNode%p_Fsingle1D(1:isizeCopy))
        case (ST_DOUBLE)
          call lalg_copyVectorDble (p_rnode%p_Ddouble1D(1:isizeCopy),&
                                    rstorageNode%p_Ddouble1D(1:isizeCopy))
        case (ST_INT)
          call lalg_copyVectorInt (p_rnode%p_Iinteger1D(1:isizeCopy),&
                                  rstorageNode%p_Iinteger1D(1:isizeCopy))
        case (ST_LOGICAL)
          ! Copy by hand
          do i = 1, isizeCopy
            rstorageNode%p_Blogical1D(i) = p_rnode%p_Blogical1D(i)
          end do
        case (ST_CHAR)
          ! Copy by hand
          do i = 1, isizeCopy
            rstorageNode%p_Schar1D(i) = p_rnode%p_Schar1D(i)
          end do
        end select
      end if

    case (2)

      ! Get the size of the old storage node.
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        Isize2Dold = shape(p_rnode%p_Fsingle2D)
      case (ST_DOUBLE)
        Isize2Dold = shape(p_rnode%p_Ddouble2D)
      case (ST_INT)
        Isize2Dold = shape(p_rnode%p_Iinteger2D)
      case (ST_LOGICAL)
        Isize2Dold = shape(p_rnode%p_Blogical2D)
      case (ST_CHAR)
        Isize2Dold = shape(p_rnode%p_Schar2D)
      end select

      ! Do we really have to change anything?
      if (isize == Isize2Dold(2)) return

      ! Allocate new memory and initialise it - if it's larger than the old
      ! memory block.
      select case (rstorageNode%idataType)
      case (ST_SINGLE)
        allocate(rstorageNode%p_Fsingle2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_SINGLE2BYTES,DP)
      case (ST_DOUBLE)
        allocate(rstorageNode%p_Ddouble2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_DOUBLE2BYTES,DP)
      case (ST_INT)
        allocate(rstorageNode%p_Iinteger2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_INT2BYTES,DP)
      case (ST_LOGICAL)
        allocate(rstorageNode%p_Blogical2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_LOGICAL2BYTES,DP)
      case (ST_CHAR)
        allocate(rstorageNode%p_Schar2D(Isize2Dold(1),isize))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_CHAR2BYTES,DP)
      case DEFAULT
        print *,'Error: unknown mem type'
        call sys_halt()
      end select

      if (isize > Isize2Dold(2)) &
        call storage_initialiseNode (rstorageNode,cinitNewBlock,&
                                     Isize2Dold(2)+1_I32)

      ! Copy old data?
      if (bcopyData) then

        ! Here it's easier than in storage_copy as we can be sure, source and
        ! destination array have the same type!
        select case (rstorageNode%idataType)
        case (ST_DOUBLE)
          ! Copy by hand
          do j=1,min(size(rstorageNode%p_Ddouble2D,2),Isize2DOld(2))
            do i=1,size(rstorageNode%p_Ddouble2D,1)
              rstorageNode%p_Ddouble2D(i,j) = p_rnode%p_Ddouble2D(i,j)
            end do
          end do

        case (ST_SINGLE)
          ! Copy by hand
          do j=1,min(size(rstorageNode%p_Fsingle2D,2),Isize2DOld(2))
            do i=1,size(rstorageNode%p_Fsingle2D,1)
              rstorageNode%p_Fsingle2D(i,j) = p_rnode%p_Fsingle2D(i,j)
            end do
          end do

        case (ST_INT)
          ! Copy by hand
          do j=1,min(size(rstorageNode%p_Iinteger2D,2),Isize2DOld(2))
            do i=1,size(rstorageNode%p_Iinteger2D,1)
              rstorageNode%p_Iinteger2D(i,j) = p_rnode%p_Iinteger2D(i,j)
            end do
          end do

        case (ST_LOGICAL)
          ! Copy by hand
          do j=1,min(size(rstorageNode%p_Blogical2D,2),Isize2DOld(2))
            do i=1,size(rstorageNode%p_Blogical2D,1)
              rstorageNode%p_Blogical2D(i,j) = p_rnode%p_Blogical2D(i,j)
            end do
          end do

        case (ST_CHAR)
          ! Copy by hand
          do j=1,min(size(rstorageNode%p_Schar2D,2),Isize2DOld(2))
            do i=1,size(rstorageNode%p_Schar2D,1)
              rstorageNode%p_Schar2D(i,j) = p_rnode%p_Schar2D(i,j)
            end do
          end do
        end select

      end if

    case DEFAULT
      print *, 'Error in storage_realloc: Handle ',ihandle, &
               ' is neither 1- nor 2- dimensional!'
      call sys_halt()

    end select

    ! Respect also the temporary memory in the total amount of memory used.
    if ((p_rheap%dtotalMem + rstorageNode%dmemBytes) .gt. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem + rstorageNode%dmemBytes

    ! Release old data
    if (associated(p_rnode%p_Fsingle1D))  deallocate(p_rnode%p_Fsingle1D)
    if (associated(p_rnode%p_Ddouble1D))  deallocate(p_rnode%p_Ddouble1D)
    if (associated(p_rnode%p_Iinteger1D)) deallocate(p_rnode%p_Iinteger1D)
    if (associated(p_rnode%p_Blogical1D)) deallocate(p_rnode%p_Blogical1D)
    if (associated(p_rnode%p_Schar1D))    deallocate(p_rnode%p_Schar1D)
    if (associated(p_rnode%p_Fsingle2D))  deallocate(p_rnode%p_Fsingle2D)
    if (associated(p_rnode%p_Ddouble2D))  deallocate(p_rnode%p_Ddouble2D)
    if (associated(p_rnode%p_Iinteger2D)) deallocate(p_rnode%p_Iinteger2D)
    if (associated(p_rnode%p_Blogical2D)) deallocate(p_rnode%p_Blogical2D)
    if (associated(p_rnode%p_Schar2D))    deallocate(p_rnode%p_Schar2D)

    ! Correct the memory statistics
    p_rheap%dtotalMem = p_rheap%dtotalMem &
                      - p_rnode%dmemBytes + rstorageNode%dmemBytes
    if (p_rheap%dtotalMem .gt. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem

    ! Replace the old node by the new one, finish
    p_rnode = rstorageNode

  end subroutine storage_realloc

!************************************************************************

!<subroutine>

  subroutine storage_reallocFixed (scall, ilbound, iubound, ihandle, &
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
  character(LEN=*), intent(IN) :: scall

  ! Requested lower bound for the memory block / the new lower bound of the last
  ! dimension in the memory block identified by ihandle
  integer(I32), intent(IN) :: ilbound

  ! Requested upper bound for the memory block / the new upper bound of the last
  ! dimension in the memory block identified by ihandle
  integer(I32), intent(IN) :: iubound

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO,
  ! ST_NEWBLOCK_NOINIT, ST_NEWBLOCK_ORDERED).
  ! Specifies how to initialise memory if isize > original array size.
  integer, intent(IN) :: cinitNewBlock

  ! OPTIONAL: Copy old data.
  ! =TRUE: Copy data of old array to the new one.
  ! =FALSE: Reallocate memory, don't copy old data.
  ! If not specified, TRUE is assumed.
  logical, intent(IN), optional :: bcopy

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(INOUT), target, optional :: rheap

  ! Handle of the memory block.
  integer, intent(INOUT) :: ihandle

!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode

    ! New storage node
    type(t_storageNode) :: rstorageNode

    ! size of the new 1-dimensional array
    integer(I32) :: isize

    ! size of the old 1-dimensional array
    integer(I32) :: isizeOld

    ! lower bound of the old 1-dimensional array
    integer(I32) :: ilboundOld

    ! upper bound of the old 1-dimensional array
    integer(I32) :: iuboundOld

    ! lower bound of the 1-dimensional array to be copied
    integer :: ilboundCopy

    ! upper bound of the 1-dimensional array to be copied
    integer :: iuboundCopy

    ! size of the old 2-dimensional array
    integer(I32), dimension(2) :: Isize2Dold

    ! lower bound of the old 2-dimensional array
    integer(I32), dimension(2) :: ilbound2Dold

    ! upper bound of the old 2-dimensional array
    integer(I32), dimension(2) :: iubound2Dold

    integer(I32) :: i,j

    logical :: bcopyData

    isize=iubound-ilbound+1
    if (isize .eq. 0) then
      ! Ok, not much to do...
      call storage_free(ihandle,rheap)
      return
    end if

    ! Get the heap to use - local or global one.

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    ! Copy old data?

    if (present(bcopy)) then
      bcopyData = bcopy
    else
      bcopyData = .true.
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Copy the data of the old storage node to rstorageNode. That way
    ! we prepare a new storage node and will replace the old.
    rstorageNode = p_rnode

    ! Are we 1D or 2D?
    select case(p_rnode%idimension)

    case (1)

      ! Get the size and bounds of the old storage node.
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        isizeOld   = size(p_rnode%p_Fsingle1D)
        ilboundOld = lbound(p_rnode%p_Fsingle1D,1)
        iuboundOld = ubound(p_rnode%p_Fsingle1D,1)
      case (ST_DOUBLE)
        isizeOld   = size(p_rnode%p_Ddouble1D)
        ilboundOld = lbound(p_rnode%p_Ddouble1D,1)
        iuboundOld = ubound(p_rnode%p_Ddouble1D,1)
      case (ST_INT)
        isizeOld   = size(p_rnode%p_Iinteger1D)
        ilboundOld = lbound(p_rnode%p_Iinteger1D,1)
        iuboundOld = ubound(p_rnode%p_Iinteger1D,1)
      case (ST_LOGICAL)
        isizeOld   = size(p_rnode%p_Blogical1D)
        ilboundOld = lbound(p_rnode%p_Blogical1D,1)
        iuboundOld = ubound(p_rnode%p_Blogical1D,1)
      case (ST_CHAR)
        isizeOld   = size(p_rnode%p_Schar1D)
        ilboundOld = lbound(p_rnode%p_Schar1D,1)
        iuboundOld = ubound(p_rnode%p_Schar1D,1)
      end select

      ! Do we really have to change anything?
      if ((ilbound == ilboundOld) .and. &
          (iubound == iuboundOld)) return

      ! Allocate new memory and initialise it - if it's larger than the old
      ! memory block.

      select case (rstorageNode%idataType)
      case (ST_SINGLE)
        allocate(rstorageNode%p_Fsingle1D(ilbound:iubound))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_SINGLE2BYTES,DP)
      case (ST_DOUBLE)
        allocate(rstorageNode%p_Ddouble1D(ilbound:iubound))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_DOUBLE2BYTES,DP)
      case (ST_INT)
        allocate(rstorageNode%p_Iinteger1D(ilbound:iubound))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_INT2BYTES,DP)
      case (ST_LOGICAL)
        allocate(rstorageNode%p_Blogical1D(ilbound:iubound))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_LOGICAL2BYTES,DP)
      case (ST_CHAR)
        allocate(rstorageNode%p_Schar1D(ilbound:iubound))
        rstorageNode%dmemBytes = real(isize,DP)*real(ST_CHAR2BYTES,DP)
      case DEFAULT
        print *,'Error: unknown mem type'
        call sys_halt()
      end select

      if (iubound > iuboundOld) &
          call storage_initialiseNode (rstorageNode,cinitNewBlock,iuboundOld+1_I32)
      if (ilbound < ilboundOld) &
          call storage_initialiseNode (rstorageNode,cinitNewBlock,ilbound,ilboundOld-1_I32)

      ! Copy old data?
      if (bcopyData) then
        ilboundCopy=max(ilbound,ilboundOld)
        iuboundCopy=min(iubound,iuboundOld)
        select case (rstorageNode%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorSngl (p_rnode%p_Fsingle1D(ilboundCopy:iuboundCopy),&
                                    rstorageNode%p_Fsingle1D(ilboundCopy:iuboundCopy))
        case (ST_DOUBLE)
          call lalg_copyVectorDble (p_rnode%p_Ddouble1D(ilboundCopy:iuboundCopy),&
                                    rstorageNode%p_Ddouble1D(ilboundCopy:iuboundCopy))
        case (ST_INT)
          call lalg_copyVectorInt (p_rnode%p_Iinteger1D(ilboundCopy:iuboundCopy),&
                                  rstorageNode%p_Iinteger1D(ilboundCopy:iuboundCopy))
        case (ST_LOGICAL)
          ! Copy by hand
          do i = ilboundCopy, iuboundCopy
            rstorageNode%p_Blogical1D(i) = p_rnode%p_Blogical1D(i)
          end do
        case (ST_CHAR)
          ! Copy by hand
          do i = ilboundCopy, iuboundCopy
            rstorageNode%p_Schar1D(i) = p_rnode%p_Schar1D(i)
          end do
        end select
      end if

    case (2)

      ! Get the size of the old storage node.
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        Isize2Dold = shape(p_rnode%p_Fsingle2D)
        ilbound2Dold = lbound(p_rnode%p_Fsingle2D)
        iubound2Dold = ubound(p_rnode%p_Fsingle2D)
      case (ST_DOUBLE)
        Isize2Dold = shape(p_rnode%p_Ddouble2D)
        ilbound2Dold = lbound(p_rnode%p_Ddouble2D)
        iubound2Dold = ubound(p_rnode%p_Ddouble2D)
      case (ST_INT)
        Isize2Dold = shape(p_rnode%p_Iinteger2D)
        ilbound2Dold = lbound(p_rnode%p_Iinteger2D)
        iubound2Dold = ubound(p_rnode%p_Iinteger2D)
      case (ST_LOGICAL)
        Isize2Dold = shape(p_rnode%p_Blogical2D)
        ilbound2Dold = lbound(p_rnode%p_Blogical2D)
        iubound2Dold = ubound(p_rnode%p_Blogical2D)
      case (ST_CHAR)
        Isize2Dold = shape(p_rnode%p_Schar2D)
        ilbound2Dold = lbound(p_rnode%p_Schar2D)
        iubound2Dold = ubound(p_rnode%p_Schar2D)
      end select

      ! Do we really have to change anything?
      if ((ilbound == ilbound2Dold(2)) .and. &
          (iubound == iubound2Dold(2))) return

      ! Allocate new memory and initialise it - if it's larger than the old
      ! memory block.
      select case (rstorageNode%idataType)
      case (ST_SINGLE)
        allocate(rstorageNode%p_Fsingle2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            &ilbound:iubound))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_SINGLE2BYTES,DP)
      case (ST_DOUBLE)
        allocate(rstorageNode%p_Ddouble2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            ilbound:iubound))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_DOUBLE2BYTES,DP)
      case (ST_INT)
        allocate(rstorageNode%p_Iinteger2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            ilbound:iubound))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_INT2BYTES,DP)
      case (ST_LOGICAL)
        allocate(rstorageNode%p_Blogical2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            ilbound:iubound))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_LOGICAL2BYTES,DP)
      case (ST_CHAR)
        allocate(rstorageNode%p_Schar2D(&
            Ilbound2Dold(1):Iubound2Dold(1),&
            ilbound:iubound))
        rstorageNode%dmemBytes = &
             real(Isize2Dold(1),DP)*real(isize,DP)*real(ST_CHAR2BYTES,DP)
      case DEFAULT
        print *,'Error: unknown mem type'
        call sys_halt()
      end select

      if (iubound > Iubound2Dold(2)) &
        call storage_initialiseNode (rstorageNode,cinitNewBlock,&
                                     Iubound2Dold(2)+1_I32)
      if (ilbound < Ilbound2Dold(2)) &
        call storage_initialiseNode (rstorageNode,cinitNewBlock,&
                                     ilbound,Ilbound2Dold(2)-1_I32)

      ! Copy old data?
      if (bcopyData) then

        ! Here it's easier than in storage_copy as we can be sure, source and
        ! destination array have the same type!
        select case (rstorageNode%idataType)
        case (ST_DOUBLE)
          ! Copy by hand
          do j=max(ilbound,Ilbound2Dold(2)),min(iubound,Iubound2Dold(2))
            do i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Ddouble2D(i,j) = p_rnode%p_Ddouble2D(i,j)
            end do
          end do

        case (ST_SINGLE)
          ! Copy by hand
          do j=max(ilbound,Ilbound2Dold(2)),min(iubound,Iubound2Dold(2))
            do i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Fsingle2D(i,j) = p_rnode%p_Fsingle2D(i,j)
            end do
          end do

        case (ST_INT)
          ! Copy by hand
          do j=max(ilbound,Ilbound2Dold(2)),min(iubound,Iubound2Dold(2))
            do i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Iinteger2D(i,j) = p_rnode%p_Iinteger2D(i,j)
            end do
          end do

        case (ST_LOGICAL)
          ! Copy by hand
          do j=max(ilbound,Ilbound2Dold(2)),min(iubound,Iubound2Dold(2))
            do i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Blogical2D(i,j) = p_rnode%p_Blogical2D(i,j)
            end do
          end do

        case (ST_CHAR)
          ! Copy by hand
          do j=max(ilbound,Ilbound2Dold(2)),min(iubound,Iubound2Dold(2))
            do i=Ilbound2DOld(1),Iubound2Dold(1)
              rstorageNode%p_Schar2D(i,j) = p_rnode%p_Schar2D(i,j)
            end do
          end do
        end select

      end if

    case DEFAULT
      print *, 'Error in storage_realloc: Handle ',ihandle, &
               ' is neither 1- nor 2- dimensional!'
      call sys_halt()

    end select

    ! Respect also the temporary memory in the total amount of memory used.
    if ((p_rheap%dtotalMem + rstorageNode%dmemBytes) .gt. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem + rstorageNode%dmemBytes

    ! Release old data
    if (associated(p_rnode%p_Fsingle1D))  deallocate(p_rnode%p_Fsingle1D)
    if (associated(p_rnode%p_Ddouble1D))  deallocate(p_rnode%p_Ddouble1D)
    if (associated(p_rnode%p_Iinteger1D)) deallocate(p_rnode%p_Iinteger1D)
    if (associated(p_rnode%p_Blogical1D)) deallocate(p_rnode%p_Blogical1D)
    if (associated(p_rnode%p_Schar1D))    deallocate(p_rnode%p_Schar1D)
    if (associated(p_rnode%p_Fsingle2D))  deallocate(p_rnode%p_Fsingle2D)
    if (associated(p_rnode%p_Ddouble2D))  deallocate(p_rnode%p_Ddouble2D)
    if (associated(p_rnode%p_Iinteger2D)) deallocate(p_rnode%p_Iinteger2D)
    if (associated(p_rnode%p_Blogical2D)) deallocate(p_rnode%p_Blogical2D)
    if (associated(p_rnode%p_Schar2D))    deallocate(p_rnode%p_Schar2D)

    ! Correct the memory statistics
    p_rheap%dtotalMem = p_rheap%dtotalMem &
                      - p_rnode%dmemBytes + rstorageNode%dmemBytes
    if (p_rheap%dtotalMem .gt. p_rheap%dtotalMemMax) &
      p_rheap%dtotalMemMax = p_rheap%dtotalMem

    ! Replace the old node by the new one, finish
    p_rnode = rstorageNode

  end subroutine storage_reallocFixed

!************************************************************************

!<function>

  function storage_isEqual (ihandle1, ihandle2, rheap1, rheap2) result (bisequal)

!<description>

  ! This function checks if the content of two different handles is equal

!</description>

!<input>

  ! The first handle
  integer, intent(IN) :: ihandle1

  ! The second handle
  integer, intent(IN) :: ihandle2

  ! OPTIONAL: first local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap1

  ! OPTIONAL: second local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(IN), target, optional :: rheap2

!</input>

!<result>

  ! .TRUE. if the content is equal.
  logical :: bisequal

!</result>

!</function>

  ! local variables

  ! Pointer to the heaps
  type(t_storageBlock), pointer :: p_rheap1,p_rheap2

  ! Identifier for data type
  integer :: idatatype1, idatatype2

  ! Identifier for data dimension
  integer :: idimension1, idimension2

  ! Identifier for size
  integer(I32) :: isize1, isize2
  integer(I32) :: Isize2D1(2), Isize2D2(2)

  ! Auxiliary arrays
  real(DP), dimension(:,:), pointer     :: p_Ddouble2D1,p_Ddouble2D2
  real(DP), dimension(:),   pointer     :: p_Ddouble1D1,p_Ddouble1D2
  real(SP), dimension(:,:), pointer     :: p_Fsingle2D1,p_Fsingle2D2
  real(SP), dimension(:),   pointer     :: p_Fsingle1D1,p_Fsingle1D2
  integer(I32), dimension(:,:), pointer :: p_Iinteger2D1,p_Iinteger2D2
  integer(I32), dimension(:),   pointer :: p_Iinteger1D1,p_Iinteger1D2
  logical, dimension(:,:), pointer      :: p_Blogical2D1,p_Blogical2D2
  logical, dimension(:),   pointer      :: p_Blogical1D1,p_Blogical1D2
  character, dimension(:,:), pointer    :: p_Schar2D1,p_Schar2D2
  character, dimension(:),   pointer    :: p_Schar1D1,p_Schar1D2
  
  ! Get the heaps to use - local or global one.

  if(present(rheap1)) then
    p_rheap1 => rheap1
  else
    p_rheap1 => rbase
  end if

  if(present(rheap2)) then
    p_rheap2 => rheap2
  else
    p_rheap2 => rbase
  end if

  ! Determine data type
  call storage_getdatatype(ihandle1, idatatype1, p_rheap1)
  call storage_getdatatype(ihandle2, idatatype2, p_rheap2)

  if (idatatype1 .ne. idatatype2) then
    bisequal = .false.
    return
  end if

  ! Determine dimension
  call storage_getdimension(ihandle1, idimension1, p_rheap1)
  call storage_getdimension(ihandle2, idimension2, p_rheap2)

  if (idimension1 .ne. idimension2) then
    bisequal = .false.
    return
  end if

  ! What dimension do we have?
  select case(idimension1)
  case (1)

    ! Determine size
    call storage_getsize1d(ihandle1, isize1, p_rheap1)
    call storage_getsize1d(ihandle2, isize2, p_rheap2)

    if (isize1 .ne. isize2) then
      bisequal = .false.
      return
    end if

    ! What data type do we have?
    select case(idatatype1)
    case (ST_INT)
      call storage_getbase_int(ihandle1, p_Iinteger1D1, p_rheap1)
      call storage_getbase_int(ihandle2, p_Iinteger1D2, p_rheap2)

      bisequal = all(p_Iinteger1D1 .eq. p_Iinteger1D2)

    case (ST_SINGLE)
      call storage_getbase_single(ihandle1, p_Fsingle1D1, p_rheap1)
      call storage_getbase_single(ihandle2, p_Fsingle1D2, p_rheap2)
      
      bisequal = all(p_Fsingle1D1 .eq. p_Fsingle1D2)

    case (ST_DOUBLE)
      call storage_getbase_double(ihandle1, p_Ddouble1D1, p_rheap1)
      call storage_getbase_double(ihandle2, p_Ddouble1D2, p_rheap2)

      bisequal = all(p_Ddouble1D1 .eq. p_Ddouble1D2)
      
    case (ST_LOGICAL)
      call storage_getbase_logical(ihandle1, p_Blogical1D1, p_rheap1)
      call storage_getbase_logical(ihandle2, p_Blogical1D2, p_rheap2)

      bisequal = all(p_Blogical1D1 .and. p_Blogical1D2)
      
    case (ST_CHAR)
      call storage_getbase_char(ihandle1, p_Schar1D1, p_rheap1)
      call storage_getbase_char(ihandle2, p_Schar1D2, p_rheap2)

      bisequal = all(p_Schar1D1 .eq. p_Schar1D2)
    end select

  case (2)

    ! Determine size
    call storage_getsize2d(ihandle1, Isize2D1, p_rheap1)
    call storage_getsize2d(ihandle2, Isize2D2, p_rheap2)

    if (any(Isize2d1 .ne. Isize2D2)) then
      bisequal = .false.
      return
    end if

    ! What data type do we have?
    select case(idatatype1)
    case (ST_INT)
      call storage_getbase_int2d(ihandle1, p_Iinteger2D1, p_rheap1)
      call storage_getbase_int2d(ihandle2, p_Iinteger2D2, p_rheap2)

      bisequal = all(p_Iinteger2D1 .eq. p_Iinteger2D2)

    case (ST_SINGLE)
      call storage_getbase_single2d(ihandle1, p_Fsingle2D1, p_rheap1)
      call storage_getbase_single2d(ihandle2, p_Fsingle2D2, p_rheap2)
      
      bisequal = all(p_Fsingle2D1 .eq. p_Fsingle2D2)

    case (ST_DOUBLE)
      call storage_getbase_double2d(ihandle1, p_Ddouble2D1, p_rheap1)
      call storage_getbase_double2d(ihandle2, p_Ddouble2D2, p_rheap2)

      bisequal = all(p_Ddouble2D1 .eq. p_Ddouble2D2)
      
    case (ST_LOGICAL)
      call storage_getbase_logical2d(ihandle1, p_Blogical2D1, p_rheap1)
      call storage_getbase_logical2d(ihandle2, p_Blogical2D2, p_rheap2)

      bisequal = all(p_Blogical2D1 .and. p_Blogical2D2)
      
    case (ST_CHAR)
      call storage_getbase_char2d(ihandle1, p_Schar2D1, p_rheap1)
      call storage_getbase_char2d(ihandle2, p_Schar2D2, p_rheap2)

      bisequal = all(p_Schar2D1 .eq. p_Schar2D2)
    end select

  case DEFAULT
    print *, 'Error in storage_realloc: Handle ',ihandle1, &
        ' is neither 1- nor 2- dimensional!'
    call sys_halt()
  end select
  end function storage_isEqual

!************************************************************************

!<function>

  function storage_isFree (ihandle, rheap) result (bisfree)

!<description>

    ! This function checks if the handle is free

!</description>

!<input>

    ! The handle
    integer, intent(IN) :: ihandle
    
    ! OPTIONAL: local heap structure to initialise. If not given, the
    ! global heap is used.
    type(t_storageBlock), intent(IN), target, optional :: rheap

!</input>

!<result>

    ! .TRUE. if the handle is free
    logical :: bisfree

!</result>

!</function>

    ! Pointer to the heaps
    type(t_storageBlock), pointer :: p_rheap
    
    ! Get the heaps to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_isFree')
      call sys_halt()
    end if

    ! Check of the handle is free
    bisfree = (p_rheap%p_Rdescriptors(ihandle)%idataType .eq. ST_NOHANDLE)
    
  end function storage_isFree

!************************************************************************

!<subroutine>

  subroutine storage_getBlock (ihandlesInUse, nhandlesTotal, ihandlesDelta,&
                               dtotalMem, nhandlesInUseMax, dtotalMemMax, rheap)

!<subroutine>

    ! This subroutine allows to get internal structures of the storage block.

!</description>

!<input>
    
    ! OPTIONAL: local heap structure to initialise. If not given, the
    ! global heap is used.
    type(t_storageBlock), intent(IN), target, optional :: rheap

!</input>

!<output>

    ! Number of handles in use
    integer, intent(OUT), optional :: ihandlesInUse

    ! Total number of handles maintained by this block
    integer, intent(OUT), optional :: nhandlesTotal

    ! Number of handles to add if there are not enough free handles
    integer, intent(OUT), optional :: ihandlesDelta

    ! Total amount of memory (in bytes) that is in use.
    real(DP), intent(OUT), optional :: dtotalMem

    ! Maximum number of handles that were in use over 
    ! the whole lifetime of this structure.
    integer, intent(OUT), optional :: nhandlesInUseMax

    ! Maximum amount of memory that was in use over the whole lifetime
    ! of this structure.
    real(DP), intent(OUT), optional :: dtotalMemMax

!</output>

!</subroutine>

    ! Pointer to the heaps
    type(t_storageBlock), pointer :: p_rheap
    
    ! Get the heaps to use - local or global one.
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (present(ihandlesInUse))    ihandlesInUse = p_rheap%ihandlesInUse
    if (present(nhandlesTotal))    nhandlesTotal = p_rheap%nhandlesTotal
    if (present(ihandlesDelta))    ihandlesDelta = p_rheap%ihandlesDelta
    if (present(dtotalMem))        dtotalMem = p_rheap%dtotalMem
    if (present(nhandlesInUseMax)) nhandlesInUseMax = p_rheap%nhandlesInUseMax
    if (present(dtotalMemMax))     dtotalMemMax = p_rheap%dtotalMemMax

  end subroutine storage_getBlock
end module storage
