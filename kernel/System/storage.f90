!##############################################################################
!# ****************************************************************************
!# <name> storage </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the FEAT implementation of the global memory
!# management. The usual memory management in Fortran 90 uses pointers
!# for referring to a memory block. As this handling is rather nasty in
!# many circumstances (e.g. it is hard to set up something like an
!# 'array of pointers'), we decided to implement our own memory management -
!# based on ALLOCATE and DEALLOCATE.
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
!# The memory management supports 1D and 2D arrays for SINGLE, DOUBLE
!# and (on some architectures) QUAD PRECISION floating point variables
!# as well as I8, I16, I32 and I64 integer variables and LOGICAL and
!# CHAR variables.
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
!#  4.) storage_new
!#      -> Allocates a new 1D, 2D or 3D array
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
!#         integers, logicals or chars
!#
!#  7.) storage_getbase_single2D,
!#      storage_getbase_double2D,
!#      storage_getbase_int2D,
!#      storage_getbase_logical2D,
!#      storage_getbase_char2D,
!#      -> Determine pointer associated to a handle for singles, doubles,
!#         integers, logicals or chars, 2D array
!#
!#  8.) storage_getbase_single3D,
!#      storage_getbase_double3D,
!#      storage_getbase_int3D,
!#      storage_getbase_logical3D,
!#      storage_getbase_char3D,
!#      -> Determine pointer associated to a handle for singles, doubles,
!#         integers, logicals or chars, 3D array
!#
!#  9.) storage_copy = storage_copy / storage_copy_explicit / storage_copy_explicit1D
!#      -> Copies the content of one array to another.
!#
!# 10.) storage_clear
!#      -> Clears an array by overwriting the entries with 0 (or .FALSE. for LOGICALs).
!#
!# 11.) storage_getsize = storage_getsize1d / storage_getsize
!#      -> Get the length of an array on the heap.
!#
!# 12.) storage_getdatatype
!#      -> Get the datatype of an array on the heap.
!#
!# 13.) storage_getdimension
!#      -> Get the dimension of an array on the heap.
!#
!# 14.) storage_realloc
!#      -> Reallocate an array (last dimension only)
!#
!# 15.) storage_initialiseBlock
!#      -> Initialise a storage block with zero (like storage_clear)
!#         or increasing number.
!#
!# 16.) storage_isEqual
!#      -> Checks if the content of two different handles is equal
!#
!# 17.) storage_createFpdbObject
!#      -> Creates an ObjectItem representing the storage management
!#
!# 18.) storage_restoreFpdbObject
!#      -> Restores the storage management from an ObjectItem
!#
!# 19.) storage_setdatatype
!#      -> Change the data type of a memory block behind a handle
!#
!# 20.) storage_getblocktype
!#      -> Get the type of the data type from a given string
!#         representation or return ST_NOHANDLE if data type is not supported
!#
!# 21.) storage_allocMemoryOnDevice
!#      -> Allocates a memory block in the device memory
!#
!# 22.) storage_deallocOnDevice
!#      -> Deallocates a memory block in the device memory
!#
!# 23.) storage_syncMemoryHostDevice
!#      -> Synchronises memory blocks between host and device memory.
!#
!# 24.) storage_clearMemoryOnDevice
!#      -> Clears a memory block in the device memory
!#
!# 25.) storage_getMemPtrOnDevice
!#      -> Get the memory address in the device memory
!# </purpose>
!##############################################################################

module storage
  
  use fsystem
  use fpersistence
  use genoutput
  use linearalgebra
  use uuid

  ! By default memory is allocated and deallocated using the Fortran
  ! intrinsics ALLOCATE/DEALLOCATE. If the Fortran compiler support
  ! ISO_C_BINDING, then any C-type malloc/free can be used to allocate
  ! and deallocate memory and associated it with the memory block
#define C_PTR_STORAGE_MALLOC 1
#define C_PTR_STORAGE_COPROC 2

  !*****************************************************************************
  ! Check for ISO_C_BINDING support
  !*****************************************************************************

#ifdef HAS_ISO_C_BINDING
  ! The external module file iso_c_binding cannot be included by the
  ! standard use statement since the configure script would generate a
  ! rule for building iso_c_binding.mod from iso_c_binfing.f90 which,
  ! of cource, does not exists. Thus, the module is hidden from configure.
#define __external_use__(module) use module
  __external_use__(iso_c_binding)

#else

#ifdef USE_C_PTR_STORAGE
  ! If the compiler does not support iso_c_binding but use of
  ! C_PTR_STORAGE has been defined then throw an error and stop
#error "Compiler does not support ISO_C_BINGING!"
#endif
 
#endif

  implicit none
  
  private

#ifndef HAS_ISO_C_BINDING
  ! If the compiler does not provice ISO_C_BINDING we need to define
  ! C_PTR and C_NULL_PTR by hand so that they can be exported below
  type C_PTR
    private
#if defined(_WIN64) || defined(_LP64) || defined(__LP64__)
    ! LP64 machine, OS X or Linux or Unix or
    ! LLP64 machine, Windows
    integer(I64) :: imemAddress = 0_I64
#elif defined(_WIN32) || defined(_ILD32)
    ! 32-bit machine, Windows or Linux or OS X or Unix
    integer(I32) :: imemAddress = 0_I32
#else
    ! Machine not detected uniquely. Since most machines are 64-bit
    ! nowadays, C_PTR is defined as 64-bit pointer
    integer(I64) :: imemAddress = 0_I64
#endif
  end type C_PTR

  ! Null pointer
  type(C_PTR), parameter :: C_NULL_PTR = C_PTR(0)
#endif
  
  public :: C_PTR

!<constants>

!<constantblock description="Storage block type identifiers">

  ! defines a non-allocated storage block handle
  integer, parameter, public :: ST_NOHANDLE = 0

  ! storage block contains single precision floats
  integer, parameter, public :: ST_SINGLE = 1

  ! storage block contains double precision floats
  integer, parameter, public :: ST_DOUBLE = 2
  
  ! storage block contains quad precision floats
  integer, parameter, public :: ST_QUAD = 3

  ! storage block contains integers (compiler-dependent size)
  integer, parameter, public :: ST_INT = 4

  ! storage block contains 8 bit integers
  integer, parameter, public :: ST_INT8 = 5

  ! storage block contains 16 bit integers
  integer, parameter, public :: ST_INT16 = 6

  ! storage block contains 32 bit integers
  integer, parameter, public :: ST_INT32 = 7
  
  ! storage block contains 64 bit integers
  integer, parameter, public :: ST_INT64 = 8

  ! storage block contains logicals
  integer, parameter, public :: ST_LOGICAL = 9

  ! storage block contains characters
  integer, parameter, public :: ST_CHAR = 10

!</constantblock>

!<constantblock description="Constants for initialisation of memory on allocation">

  ! init new storage block with zeros (or .FALSE. for logicals)
  integer, parameter, public :: ST_NEWBLOCK_ZERO = 0

  ! no init new storage block
  integer, parameter, public :: ST_NEWBLOCK_NOINIT = 1

  ! init new storage block with 1,2,3,...,n
  ! Note: This initialization constant must NOT be used for initialisation of logicals
  !       or characters!!!
  integer, parameter, public :: ST_NEWBLOCK_ORDERED = 2

!</constantblock>

!<constantblock description="Constants for calculating memory">

  ! How many bytes has a single precision real?
  integer(I64), public :: ST_SINGLE_BYTES = int(SP,I64)

  ! How many bytes has a double precision real?
  integer(I64), public :: ST_DOUBLE_BYTES = int(DP,I64)
  
  ! How many bytes has a quad precision real?
  integer(I64), public :: ST_QUAD_BYTES = int(QP,I64)

  ! How many bytes has an integer?
  integer(I64), public :: ST_INT_BYTES = int((bit_size(1) / 8),I64)

  ! How many bytes has a 8 bit integer? 1, of course
  integer(I64), public :: ST_INT8_BYTES = int(I8,I64)

  ! How many bytes has a 16 bit integer? 2, of course
  integer(I64), public :: ST_INT16_BYTES = int(I16,I64)

  ! How many bytes has a 32 bit integer? 4, of course
  integer(I64), public :: ST_INT32_BYTES = int(I32,I64)
  
  ! How many bytes has a 64 bit integer? 8, of course
  integer(I64), public :: ST_INT64_BYTES = int(I64,I64)

  ! How many bytes has a logical?
  ! The FORTRAN standard requires logical variables to be the same size
  ! as INTEGER variables. Thus, their size may change if the size of
  ! default INTEGER variables changes, e.g. by specifying compiled switches.
  integer(I64), public :: ST_LOGICAL_BYTES = int((bit_size(1) / 8),I64)

  ! How many bytes has a character?
  ! Note: We are not 100% sure, but this may differ on other architectures... O_o
  integer(I64), public :: ST_CHAR_BYTES = 1_I64
  
!</constantblock>

!<constantblock description="Constants for synchronising memory">

  ! copy storage block from host memory to device memory (overwrite)
  integer, parameter, public :: ST_SYNCBLOCK_COPY_H2D       = 1
  
  ! copy storage block from device memory to host memory (overwrite)
  integer, parameter, public :: ST_SYNCBLOCK_COPY_D2H       = 2

  ! copy storage block from host memory to device memory (accumulate)
  integer, parameter, public :: ST_SYNCBLOCK_ACCUMULATE_H2D = 3

  ! copy storage block from device memory to host memory (accumulate)
  integer, parameter, public :: ST_SYNCBLOCK_ACCUMULATE_D2H = 4

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
    ! ST_DOUBLE, ST_INT, ST_LOGICAL, ST_CHAR)
    integer :: idataType = ST_NOHANDLE

    ! Dimension associated to the handle
    ! (0=not assigned, 1=1D, 2=2D array, 3=3D array)
    integer :: idimension = 0

    ! The name of the array that is associated to that handle
    character(LEN=SYS_NAMELEN) :: sname = ''

    ! Amount of memory (in bytes) associated to this block.
    integer(I64) :: imemBytes = 0_I64

    ! Pointer to 1D real array or NULL() if not assigned
    real(SP), dimension(:), pointer       :: p_Fsingle1D    => null()

    ! Pointer to 1D double precision array or NULL() if not assigned
    real(DP), dimension(:), pointer       :: p_Ddouble1D    => null()

    ! Pointer to 1D quad precision array or NULL() if not assigned
    real(QP), dimension(:), pointer       :: p_Qquad1D      => null()

    ! Pointer to 1D integer array or NULL() if not assigned
    integer, dimension(:), pointer        :: p_Iinteger1D   => null()
    
    ! Pointer to 1D integer(I8) array or NULL() if not assigned
    integer(I8), dimension(:), pointer    :: p_Iint8_1D     => null()
    
    ! Pointer to 1D integer(I12) array or NULL() if not assigned
    integer(I16), dimension(:), pointer   :: p_Iint16_1D    => null()
    
    ! Pointer to 1D integer(I32) array or NULL() if not assigned
    integer(I32), dimension(:), pointer   :: p_Iint32_1D    => null()

    ! Pointer to 1D integer(I64) array or NULL() if not assigned
    integer(I64), dimension(:), pointer   :: p_Iint64_1D    => null()

    ! Pointer to 1D logical array or NULL() if not assigned
    logical, dimension(:), pointer        :: p_Blogical1D   => null()

    ! Pointer to 1D character array or NULL() if not assigned
    character, dimension(:), pointer      :: p_Schar1D      => null()

    ! Pointer to 2D real array or NULL() if not assigned
    real(SP), dimension(:,:), pointer     :: p_Fsingle2D    => null()

    ! Pointer to 2D double precision array or NULL() if not assigned
    real(DP), dimension(:,:), pointer     :: p_Ddouble2D    => null()

    ! Pointer to 2D quad precision array or NULL() if not assigned
    real(QP), dimension(:,:), pointer     :: p_Qquad2D      => null()

    ! Pointer to 2D integer array or NULL() if not assigned
    integer, dimension(:,:), pointer      :: p_Iinteger2D   => null()
    
    ! Pointer to 2D integer(I8) array or NULL() if not assigned
    integer(I8), dimension(:,:), pointer  :: p_Iint8_2D     => null()
    
    ! Pointer to 2D integer(I16) array or NULL() if not assigned
    integer(I16), dimension(:,:), pointer :: p_Iint16_2D    => null()
    
    ! Pointer to 2D integer(I32) array or NULL() if not assigned
    integer(I32), dimension(:,:), pointer :: p_Iint32_2D    => null()

    ! Pointer to 2D integer(I64) array or NULL() if not assigned
    integer(I64), dimension(:,:), pointer :: p_Iint64_2D    => null()

    ! Pointer to 2D logical array or NULL() if not assigned
    logical, dimension(:,:), pointer      :: p_Blogical2D   => null()

    ! Pointer to 2D character array or NULL() if not assigned
    character, dimension(:,:), pointer    :: p_Schar2D      => null()

    ! Pointer to 3D real array or NULL() if not assigned
    real(SP), dimension(:,:,:), pointer   :: p_Fsingle3D    => null()

    ! Pointer to 3D double precision array or NULL() if not assigned
    real(DP), dimension(:,:,:), pointer   :: p_Ddouble3D    => null()

    ! Pointer to 3D quad precision array or NULL() if not assigned
    real(QP), dimension(:,:,:), pointer   :: p_Qquad3D      => null()

    ! Pointer to 3D integer array or NULL() if not assigned
    integer, dimension(:,:,:), pointer    :: p_Iinteger3D   => null()
    
    ! Pointer to 3D integer(I8) array or NULL() if not assigned
    integer(I8), dimension(:,:,:), pointer :: p_Iint8_3D    => null()
    
    ! Pointer to 3D integer(I12) array or NULL() if not assigned
    integer(I16), dimension(:,:,:), pointer :: p_Iint16_3D  => null()
    
    ! Pointer to 3D integer(I32) array or NULL() if not assigned
    integer(I32), dimension(:,:,:), pointer :: p_Iint32_3D  => null()

    ! Pointer to 3D integer(I64) array or NULL() if not assigned
    integer(I64), dimension(:,:,:), pointer :: p_Iint64_3D  => null()

    ! Pointer to 3D logical array or NULL() if not assigned
    logical, dimension(:,:,:), pointer    :: p_Blogical3D   => null()

    ! Pointer to 3D character array or NULL() if not assigned
    character, dimension(:,:,:), pointer  :: p_Schar3D      => null()
    
#ifdef USE_C_PTR_STORAGE
    ! If the code is compiled with coprocessor support and
    ! iso_c_binding enabled, then memory allocation on host and device
    ! is done by the coproc-library which supports page-locked memory.
    
    ! Pointer of type C_PTR to memory block on host
    type(C_PTR) :: chostMemPtr = C_NULL_PTR
#endif
    
#ifdef ENABLE_COPROCESSOR_SUPPORT
    ! Starting memory address on coprocessor device (=0 if not assigned)
    type(C_PTR) :: cdeviceMemPtr = C_NULL_PTR
#endif

  end type t_storageNode
  
!</typeblock>

!<typeblock>

  ! This block represents a heap that maintains single, double precision
  ! and integer data. It contains a list of t_storageNode elements for all
  ! the handles.
  ! There is one global object of this type for the global storage management,
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

    ! Total amount of memory (in bytes) that is in use.
    integer(I64) :: itotalMem = 0_I64

    ! Maximum number of handles that were in use over the whole lifetime
    ! of this structure.
    integer :: nhandlesInUseMax = 0

    ! Maximum amount of memory that was in use over the whole lifetime
    ! of this structure.
    integer(I64) :: itotalMemMax = 0_I64

  end type t_storageBlock
  
  public :: t_storageBlock

!</typeblock>

!</types>

!<publicvars>

  ! Global memory management structure
  type(t_storageBlock), private, save, target :: rbase

!</publicvars>

  interface storage_new
    module procedure storage_newDefault
    module procedure storage_newDefaultFixed
    module procedure storage_new1D
    module procedure storage_new1DFixed
    module procedure storage_newIndirect
  end interface storage_new

  public :: storage_new
  
  interface storage_getbase_int
    module procedure storage_getbase_intDefault
    module procedure storage_getbase_intUBnd
    module procedure storage_getbase_intLUBnd
  end interface
  
  public :: storage_getbase_int
  
  interface storage_getbase_int8
    module procedure storage_getbase_int8Default
    module procedure storage_getbase_int8UBnd
    module procedure storage_getbase_int8LUBnd
  end interface

  public :: storage_getbase_int8
  
  interface storage_getbase_int16
    module procedure storage_getbase_int16Default
    module procedure storage_getbase_int16UBnd
    module procedure storage_getbase_int16LUBnd
  end interface
  
  public :: storage_getbase_int16
  
  interface storage_getbase_int32
    module procedure storage_getbase_int32Default
    module procedure storage_getbase_int32UBnd
    module procedure storage_getbase_int32LUBnd
  end interface
  
  public :: storage_getbase_int32

  interface storage_getbase_int64
    module procedure storage_getbase_int64Default
    module procedure storage_getbase_int64UBnd
    module procedure storage_getbase_int64LUBnd
  end interface

  public :: storage_getbase_int64

  interface storage_getbase_single
    module procedure storage_getbase_singleDefault
    module procedure storage_getbase_singleUBnd
    module procedure storage_getbase_singleLUBnd
  end interface
  
  public :: storage_getbase_single

  interface storage_getbase_double
    module procedure storage_getbase_doubleDefault
    module procedure storage_getbase_doubleUBnd
    module procedure storage_getbase_doubleLUBnd
  end interface

  public :: storage_getbase_double

  interface storage_getbase_quad
    module procedure storage_getbase_quadDefault
    module procedure storage_getbase_quadUBnd
    module procedure storage_getbase_quadLUBnd
  end interface
  
  public :: storage_getbase_quad

  interface storage_getbase_logical
    module procedure storage_getbase_logicalDefault
    module procedure storage_getbase_logicalUBnd
    module procedure storage_getbase_logicalLUBnd
  end interface

  public :: storage_getbase_logical

  interface storage_getbase_char
    module procedure storage_getbase_charDefault
    module procedure storage_getbase_charUBnd
    module procedure storage_getbase_charLUBnd
  end interface

  public :: storage_getbase_char

  interface storage_getbase_int2D
    module procedure storage_getbase_int2DDef
    module procedure storage_getbase_int2DUBnd
    module procedure storage_getbase_int2DLUBnd
  end interface

  public :: storage_getbase_int2D
  
  interface storage_getbase_int8_2D
    module procedure storage_getbase_int8_2DDef
    module procedure storage_getbase_int8_2DUBnd
    module procedure storage_getbase_int8_2DLUBnd
  end interface
  
  public :: storage_getbase_int8_2D
  
  interface storage_getbase_int16_2D
    module procedure storage_getbase_int16_2DDef
    module procedure storage_getbase_int16_2DUBnd
    module procedure storage_getbase_int16_2DLUBnd
  end interface
  
  public :: storage_getbase_int16_2D
  
  interface storage_getbase_int32_2D
    module procedure storage_getbase_int32_2DDef
    module procedure storage_getbase_int32_2DUBnd
    module procedure storage_getbase_int32_2DLUBnd
  end interface
  
  public :: storage_getbase_int32_2D
  
  interface storage_getbase_int64_2D
    module procedure storage_getbase_int64_2DDef
    module procedure storage_getbase_int64_2DUBnd
    module procedure storage_getbase_int64_2DLUBnd
  end interface
  
  public :: storage_getbase_int64_2D

  interface storage_getbase_single2D
    module procedure storage_getbase_single2DDef
    module procedure storage_getbase_single2DUBnd
    module procedure storage_getbase_single2DLUBnd
  end interface

  public :: storage_getbase_single2D

  interface storage_getbase_double2D
    module procedure storage_getbase_double2DDef
    module procedure storage_getbase_double2DUBnd
    module procedure storage_getbase_double2DLUBnd
  end interface

  public :: storage_getbase_double2D

  interface storage_getbase_quad2D
    module procedure storage_getbase_quad2DDef
    module procedure storage_getbase_quad2DUBnd
    module procedure storage_getbase_quad2DLUBnd
  end interface

  public :: storage_getbase_quad2D

  interface storage_getbase_logical2D
    module procedure storage_getbase_logical2DDef
    module procedure storage_getbase_logical2DUBnd
    module procedure storage_getbase_logical2DLUBnd
  end interface

  public :: storage_getbase_logical2D

  interface storage_getbase_char2D
    module procedure storage_getbase_char2DDef
    module procedure storage_getbase_char2DUBnd
    module procedure storage_getbase_char2DLUBnd
  end interface

  public :: storage_getbase_char2D

  interface storage_getbase_int3D
    module procedure storage_getbase_int3DDef
    module procedure storage_getbase_int3DUBnd
    module procedure storage_getbase_int3DLUBnd
  end interface

  public :: storage_getbase_int3D
  
  interface storage_getbase_int8_3D
    module procedure storage_getbase_int8_3DDef
    module procedure storage_getbase_int8_3DUBnd
    module procedure storage_getbase_int8_3DLUBnd
  end interface
  
  public :: storage_getbase_int8_3D
  
  interface storage_getbase_int16_3D
    module procedure storage_getbase_int16_3DDef
    module procedure storage_getbase_int16_3DUBnd
    module procedure storage_getbase_int16_3DLUBnd
  end interface
  
  public :: storage_getbase_int16_3D
  
  interface storage_getbase_int32_3D
    module procedure storage_getbase_int32_3DDef
    module procedure storage_getbase_int32_3DUBnd
    module procedure storage_getbase_int32_3DLUBnd
  end interface
  
  public :: storage_getbase_int32_3D
  
  interface storage_getbase_int64_3D
    module procedure storage_getbase_int64_3DDef
    module procedure storage_getbase_int64_3DUBnd
    module procedure storage_getbase_int64_3DLUBnd
  end interface
  
  public :: storage_getbase_int64_3D

  interface storage_getbase_single3D
    module procedure storage_getbase_single3DDef
    module procedure storage_getbase_single3DUBnd
    module procedure storage_getbase_single3DLUBnd
  end interface

  public :: storage_getbase_single3D

  interface storage_getbase_double3D
    module procedure storage_getbase_double3DDef
    module procedure storage_getbase_double3DUBnd
    module procedure storage_getbase_double3DLUBnd
  end interface

  public :: storage_getbase_double3D

  interface storage_getbase_quad3D
    module procedure storage_getbase_quad3DDef
    module procedure storage_getbase_quad3DUBnd
    module procedure storage_getbase_quad3DLUBnd
  end interface

  public :: storage_getbase_quad3D

  interface storage_getbase_logical3D
    module procedure storage_getbase_logical3DDef
    module procedure storage_getbase_logical3DUBnd
    module procedure storage_getbase_logical3DLUBnd
  end interface

  public :: storage_getbase_logical3D

  interface storage_getbase_char3D
    module procedure storage_getbase_char3DDef
    module procedure storage_getbase_char3DUBnd
    module procedure storage_getbase_char3DLUBnd
  end interface

  public :: storage_getbase_char3D

  interface storage_getbase
    module procedure storage_getbase_int8Default
    module procedure storage_getbase_int8UBnd
    module procedure storage_getbase_int8LUBnd

    module procedure storage_getbase_int16Default
    module procedure storage_getbase_int16UBnd
    module procedure storage_getbase_int16LUBnd

    module procedure storage_getbase_int32Default
    module procedure storage_getbase_int32UBnd
    module procedure storage_getbase_int32LUBnd

    module procedure storage_getbase_int64Default
    module procedure storage_getbase_int64UBnd
    module procedure storage_getbase_int64LUBnd

    module procedure storage_getbase_singleDefault
    module procedure storage_getbase_singleUBnd
    module procedure storage_getbase_singleLUBnd

    module procedure storage_getbase_doubleDefault
    module procedure storage_getbase_doubleUBnd
    module procedure storage_getbase_doubleLUBnd

#ifdef ENABLE_QUADPREC
    module procedure storage_getbase_quadDefault
    module procedure storage_getbase_quadUBnd
    module procedure storage_getbase_quadLUBnd
#endif

    module procedure storage_getbase_logicalDefault
    module procedure storage_getbase_logicalUBnd
    module procedure storage_getbase_logicalLUBnd

    module procedure storage_getbase_charDefault
    module procedure storage_getbase_charUBnd
    module procedure storage_getbase_charLUBnd

    module procedure storage_getbase_int8_2DDef
    module procedure storage_getbase_int8_2DUBnd
    module procedure storage_getbase_int8_2DLUBnd

    module procedure storage_getbase_int16_2DDef
    module procedure storage_getbase_int16_2DUBnd
    module procedure storage_getbase_int16_2DLUBnd

    module procedure storage_getbase_int32_2DDef
    module procedure storage_getbase_int32_2DUBnd
    module procedure storage_getbase_int32_2DLUBnd

    module procedure storage_getbase_int64_2DDef
    module procedure storage_getbase_int64_2DUBnd
    module procedure storage_getbase_int64_2DLUBnd

    module procedure storage_getbase_single2DDef
    module procedure storage_getbase_single2DUBnd
    module procedure storage_getbase_single2DLUBnd

    module procedure storage_getbase_double2DDef
    module procedure storage_getbase_double2DUBnd
    module procedure storage_getbase_double2DLUBnd

#ifdef ENABLE_QUADPREC
    module procedure storage_getbase_quad2DDef
    module procedure storage_getbase_quad2DUBnd
    module procedure storage_getbase_quad2DLUBnd
#endif

    module procedure storage_getbase_logical2DDef
    module procedure storage_getbase_logical2DUBnd
    module procedure storage_getbase_logical2DLUBnd

    module procedure storage_getbase_char2DDef
    module procedure storage_getbase_char2DUBnd
    module procedure storage_getbase_char2DLUBnd

    module procedure storage_getbase_int8_3DDef
    module procedure storage_getbase_int8_3DUBnd
    module procedure storage_getbase_int8_3DLUBnd

    module procedure storage_getbase_int16_3DDef
    module procedure storage_getbase_int16_3DUBnd
    module procedure storage_getbase_int16_3DLUBnd

    module procedure storage_getbase_int32_3DDef
    module procedure storage_getbase_int32_3DUBnd
    module procedure storage_getbase_int32_3DLUBnd

    module procedure storage_getbase_int64_3DDef
    module procedure storage_getbase_int64_3DUBnd
    module procedure storage_getbase_int64_3DLUBnd

    module procedure storage_getbase_single3DDef
    module procedure storage_getbase_single3DUBnd
    module procedure storage_getbase_single3DLUBnd

    module procedure storage_getbase_double3DDef
    module procedure storage_getbase_double3DUBnd
    module procedure storage_getbase_double3DLUBnd

#ifdef ENABLE_QUADPREC
    module procedure storage_getbase_quad3DDef
    module procedure storage_getbase_quad3DUBnd
    module procedure storage_getbase_quad3DLUBnd
#endif

    module procedure storage_getbase_logical3DDef
    module procedure storage_getbase_logical3DUBnd
    module procedure storage_getbase_logical3DLUBnd

    module procedure storage_getbase_char3DDef
    module procedure storage_getbase_char3DUBnd
    module procedure storage_getbase_char3DLUBnd
  end interface

  public :: storage_getbase

  interface storage_getsize
    module procedure storage_getsizeDefault
    module procedure storage_getsize1D
  end interface
  
  public :: storage_getsize

  interface storage_copy
    module procedure storage_copyDefault
    module procedure storage_copy_explicit
    module procedure storage_copy_explicit1D
  end interface

  public :: storage_copy
  
  public :: storage_realloc
  public :: storage_init
  public :: storage_done
  public :: storage_free
  public :: storage_info
  public :: storage_clear
  public :: storage_getdatatype
  public :: storage_getdimension
  public :: storage_initialiseBlock
  public :: storage_isEqual
  public :: storage_createFpdbObject
  public :: storage_restoreFpdbObject
  public :: storage_setdatatype
  public :: storage_getblocktype
  public :: storage_allocMemoryOnDevice
  public :: storage_deallocMemoryOnDevice
  public :: storage_syncMemoryHostDevice
  public :: storage_clearMemoryOnDevice
  public :: storage_getMemPtrOnDevice
  public :: storage_isAssociated

  !************************************************************************

#ifdef USE_C_PTR_STORAGE

#if USE_C_PTR_STORAGE == C_PTR_STORAGE_MALLOC

  interface
    ! Define iso_c_binding Fortran: c_allocate -> C: storageHostAlloc
    integer(C_INT) function c_allocate(buffer, size) bind(C,name="storage_malloc")
      __external_use__(iso_c_binding)
      implicit none
      type (C_PTR) :: buffer
      integer (C_SIZE_T), value :: size
    end function c_allocate
  end interface
  
  interface
    ! Define iso_c_binding Fortran: c_deallocate -> C: storageFreeHost
    integer(C_INT) function c_deallocate(buffer) bind(C,name="storage_free")
      __external_use__(iso_c_binding)
      implicit none
      type (C_PTR) :: buffer
    end function c_deallocate
  end interface

#elif USE_C_PTR_STORAGE == C_PTR_STORAGE_COPROC

  interface
    ! Define iso_c_binding Fortran: c_allocate -> C: coproc_newMemoryOnHost
    integer(C_INT) function c_allocate(buffer, size) bind(C,name="coproc_malloc")
      __external_use__(iso_c_binding)
      implicit none
      type (C_PTR) :: buffer
      integer (C_SIZE_T), value :: size
    end function c_allocate
  end interface
  
  interface
    ! Define iso_c_binding Fortran: c_deallocate -> C: coproc_FreeMemoryOnHost
    integer(C_INT) function c_deallocate(buffer) bind(C,name="coproc_free")
      __external_use__(iso_c_binding)
      implicit none
      type (C_PTR) :: buffer
    end function c_deallocate
  end interface

#else
#error "Unsupported type of USE_C_PTR_STORAGE!"
#endif
#endif

contains

!************************************************************************

!<subroutine>

  subroutine storage_init (ihandleCount, ihandlesDelta, rheap)

!<description>

  ! This routine initialises the storage management.
  ! ihandleCount is the initial number of handles maintained by the
  ! storage routines. If there are not enough free handles, the number
  ! of handles are increased by ihandlesDelta (which is initially set
  ! to 1/2*ihandleCount if not given).
  ! rheap allows to specify a 'local' heap structure to initialise.
  ! If not given, the global memory management is initialised.

!</description>

!<input>

  ! Initial number of handles maintained by the storage routines.
  integer, intent(in) :: ihandleCount

  ! OPTIONAL: Number of handles to increase the memory block by, if there are
  ! not enough handles available. Standard setting is 1/2*ihandleCount.
  integer, intent(in), optional :: ihandlesDelta

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is initialised.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!</subroutine>

  ! local variables

  ! the real 'handle-delta'
  integer :: ihandles, ihDelta

  ! Pointer to the heap to initialise
  type(t_storageBlock), pointer :: p_rheap
  
  integer :: i
  integer(I64) :: isize
    
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

#ifdef ENABLE_COPROCESSOR_SUPPORT
    ! Check consistency between host and device storage managers
    call coproc_getSizeOf(ST_SINGLE, isize)
    if (isize .ne. ST_SINGLE_BYTES) then
      call output_line ('Number if bytes for ST_SINGLE mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if

    call coproc_getSizeOf(ST_DOUBLE, isize)
    if (isize .ne. ST_DOUBLE_BYTES) then
      call output_line ('Number if bytes for ST_DOUBLE mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if

    call coproc_getSizeOf(ST_QUAD, isize)
    if (isize .ne. ST_QUAD_BYTES) then
      call output_line ('Number if bytes for ST_QUAD mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if

    call coproc_getSizeOf(ST_INT, isize)
    if (isize .ne. ST_INT_BYTES) then
      call output_line ('Number if bytes for ST_INT mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if

    call coproc_getSizeOf(ST_INT8, isize)
    if (isize .ne. ST_INT8_BYTES) then
      call output_line ('Number if bytes for ST_INT8 mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if

    call coproc_getSizeOf(ST_INT16, isize)
    if (isize .ne. ST_INT16_BYTES) then
      call output_line ('Number if bytes for ST_INT16 mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if

    call coproc_getSizeOf(ST_INT32, isize)
    if (isize .ne. ST_INT32_BYTES) then
      call output_line ('Number if bytes for ST_INT32 mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if

    call coproc_getSizeOf(ST_INT64, isize)
    if (isize .ne. ST_INT64_BYTES) then
      call output_line ('Number if bytes for ST_INT64 mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if

    call coproc_getSizeOf(ST_LOGICAL, isize)
    if (isize .ne. ST_LOGICAL_BYTES) then
      call output_line ('Number if bytes for ST_LOGICAL mismatch!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_init')
      call sys_halt()
    end if
#endif
    
  end subroutine storage_init

!************************************************************************

!<subroutine>

  subroutine storage_done (rheap)

!<description>
  ! This routine cleans up the storage management. All data on the
  ! heap is released from memory.
!</description>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is cleaned up.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</inputoutput>

!</subroutine>

  ! Pointer to the heap to initialise
  type(t_storageBlock), pointer :: p_rheap
  integer :: i,ihandle
    
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    ! Delete all data from the heap
    if (associated(p_rheap%p_Rdescriptors)) then
      do i = 1,size(p_rheap%p_Rdescriptors)
        ! Do not pass i as handle as storage_free will set the handle
        ! passed to it to 0!
        ihandle = i
        if (p_rheap%p_Rdescriptors(i)%idataType .ne. ST_NOHANDLE) &
            call storage_free(ihandle,rheap)
      end do
      
      ! Release the descriptors
      deallocate(p_rheap%p_IfreeHandles)
    end if

    ! Release the descriptors
    if (associated(p_rheap%p_Rdescriptors)) then
      deallocate(p_rheap%p_Rdescriptors)
    end if
    
    ! Clean up the memory management block
    p_rheap%nhandlesTotal = 0
    p_rheap%ihandlesDelta = 0
    p_rheap%p_inextFreeHandle = 0
    p_rheap%p_ilastFreeHandle = 0
    p_rheap%ihandlesInUse = 0
        
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
  type(t_storageBlock), intent(inout) :: rheap

!</inputoutput>

!</function>

  ! local variables
  type(t_storageNode), dimension(:), pointer :: p_Rdescriptors => null()
  integer, dimension(:), pointer :: p_IfreeHandles => null()
  integer :: i
    
    ! OpenMP-Extension: It is important that no two threads generate a
    ! new handle at the same time, thus the entire routine is 'critical'

    !$omp critical(storage_global_heap_modify)

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
    
    !$omp end critical(storage_global_heap_modify)

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
  integer, intent(in) :: ihandle
!</input>

!<inputoutput>
  ! The heap structure where to release the handle from.
  type(t_storageBlock), intent(inout) :: rheap
!</inputoutput>

!</subroutine>

  type(t_storageNode), pointer :: p_rnode

    ! OpenMP-Extension: It is important that no two threads release
    ! handle at the same time, thus the entire routine is 'critical'

    !$omp critical(storage_global_heap_modify)

    ! Where is the descriptor of the handle?
    p_rnode => rheap%p_Rdescriptors(ihandle)
    
    ! Subtract the memory amount from the statistics
    rheap%itotalMem = rheap%itotalMem - p_rnode%imemBytes
    
    ! Clear the descriptor structure
    p_rnode%idataType  = ST_NOHANDLE
    p_rnode%idimension = 0
    p_rnode%imemBytes  = 0_I64
    nullify(p_rnode%p_Fsingle1D)
    nullify(p_rnode%p_Ddouble1D)
    nullify(p_rnode%p_Qquad1D)
    nullify(p_rnode%p_Iinteger1D)
    nullify(p_rnode%p_Iint8_1D)
    nullify(p_rnode%p_Iint16_1D)
    nullify(p_rnode%p_Iint32_1D)
    nullify(p_rnode%p_Iint64_1D)
    nullify(p_rnode%p_Blogical1D)
    nullify(p_rnode%p_Schar1D)

    nullify(p_rnode%p_Fsingle2D)
    nullify(p_rnode%p_Ddouble2D)
    nullify(p_rnode%p_Qquad2D)
    nullify(p_rnode%p_Iinteger2D)
    nullify(p_rnode%p_Iint8_2D)
    nullify(p_rnode%p_Iint16_2D)
    nullify(p_rnode%p_Iint32_2D)
    nullify(p_rnode%p_Iint64_2D)
    nullify(p_rnode%p_Blogical2D)
    nullify(p_rnode%p_Schar2D)

    nullify(p_rnode%p_Fsingle3D)
    nullify(p_rnode%p_Ddouble3D)
    nullify(p_rnode%p_Qquad3D)
    nullify(p_rnode%p_Iinteger3D)
    nullify(p_rnode%p_Iint8_3D)
    nullify(p_rnode%p_Iint16_3D)
    nullify(p_rnode%p_Iint32_3D)
    nullify(p_rnode%p_Iint64_3D)
    nullify(p_rnode%p_Blogical3D)
    nullify(p_rnode%p_Schar3D)
    
    ! Handle ihandle is available now - put it to the list of available handles.
    rheap%p_ilastFreeHandle = mod(rheap%p_ilastFreeHandle,rheap%nhandlesTotal) + 1
    rheap%p_IfreeHandles (rheap%p_ilastFreeHandle) = ihandle
    
    rheap%ihandlesInUse = rheap%ihandlesInUse - 1

    !$omp end critical(storage_global_heap_modify)

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
  type(t_storageNode), intent(inout) :: rstorageNode

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO, ST_NEWBLOCK_NOINIT,
  ! ST_NEWBLOCK_ORDERED). Specifies how to initialise the data block associated
  ! to ihandle.
  integer, intent(in) :: cinitNewBlock

  ! Start index from where to initialise; should usually be =1.
  ! For multidimensional arrays, this specifies the start index of the
  ! last dimension.
  integer, intent(in) :: istartIndex

  ! OPTIONAL: Stop index up to which to initialise; should usually be =SIZE(*).
  ! For multidimensional arrays, this specifies the stop index of the
  ! last dimension.
  integer, intent(in), optional :: istopIndex
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
          case (ST_QUAD)
            rstorageNode%p_Qquad1D(istartIndex:istopIndex) = 0.0_QP
          case (ST_INT)
            rstorageNode%p_Iinteger1D(istartIndex:istopIndex) = 0
          case (ST_INT8)
            rstorageNode%p_Iint8_1D(istartIndex:istopIndex) = 0_I8
          case (ST_INT16)
            rstorageNode%p_Iint16_1D(istartIndex:istopIndex) = 0_I16
          case (ST_INT32)
            rstorageNode%p_Iint32_1D(istartIndex:istopIndex) = 0_I32
          case (ST_INT64)
            rstorageNode%p_Iint64_1D(istartIndex:istopIndex) = 0_I64
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
          case (ST_QUAD)
            rstorageNode%p_Qquad1D(istartIndex:) = 0.0_QP
          case (ST_INT)
            rstorageNode%p_Iinteger1D(istartIndex:) = 0
          case (ST_INT8)
            rstorageNode%p_Iint8_1D(istartIndex:) = 0_I8
          case (ST_INT16)
            rstorageNode%p_Iint16_1D(istartIndex:) = 0_I16
          case (ST_INT32)
            rstorageNode%p_Iint32_1D(istartIndex:) = 0_I32
          case (ST_INT64)
            rstorageNode%p_Iint64_1D(istartIndex:) = 0_I64
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
          case (ST_QUAD)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Qquad1D(iorder) = real(iorder,QP)
            end do
          case (ST_INT)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Iinteger1D(iorder) = iorder
            end do
          case (ST_INT8)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Iint8_1D(iorder) = int(iorder,I8)
            end do
          case (ST_INT16)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Iint16_1D(iorder) = int(iorder,I16)
            end do
          case (ST_INT32)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Iint32_1D(iorder) = int(iorder,I32)
            end do
          case (ST_INT64)
            do iorder=istartIndex,istopIndex
              rstorageNode%p_Iint64_1D(iorder) = int(iorder,I64)
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
          case (ST_QUAD)
            do iorder=istartIndex,ubound(rstorageNode%p_Qquad1D,1)
              rstorageNode%p_Qquad1D(iorder) = real(iorder,QP)
            end do
          case (ST_INT)
            do iorder=istartIndex,ubound(rstorageNode%p_Iinteger1D,1)
              rstorageNode%p_Iinteger1D(iorder) = iorder
            end do
          case (ST_INT8)
            do iorder=istartIndex,ubound(rstorageNode%p_Iint8_1D,1)
              rstorageNode%p_Iint8_1D(iorder) = int(iorder,I8)
            end do
          case (ST_INT16)
            do iorder=istartIndex,ubound(rstorageNode%p_Iint16_1D,1)
              rstorageNode%p_Iint16_1D(iorder) = int(iorder,I16)
            end do
          case (ST_INT32)
            do iorder=istartIndex,ubound(rstorageNode%p_Iint32_1D,1)
              rstorageNode%p_Iint32_1D(iorder) = int(iorder,I32)
            end do
          case (ST_INT64)
            do iorder=istartIndex,ubound(rstorageNode%p_Iint64_1D,1)
              rstorageNode%p_Iint64_1D(iorder) = int(iorder,I64)
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
          case (ST_QUAD)
            rstorageNode%p_Qquad2D(:,istartIndex:istopIndex) = 0.0_QP
          case (ST_INT)
            rstorageNode%p_Iinteger2D(:,istartIndex:istopIndex) = 0
          case (ST_INT8)
            rstorageNode%p_Iint8_2D(:,istartIndex:istopIndex) = 0_I8
          case (ST_INT16)
            rstorageNode%p_Iint16_2D(:,istartIndex:istopIndex) = 0_I16
          case (ST_INT32)
            rstorageNode%p_Iint32_2D(:,istartIndex:istopIndex) = 0_I32
          case (ST_INT64)
            rstorageNode%p_Iint64_2D(:,istartIndex:istopIndex) = 0_I64
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
          case (ST_QUAD)
            rstorageNode%p_Qquad2D(:,istartIndex:) = 0.0_QP
          case (ST_INT)
            rstorageNode%p_Iinteger2D(:,istartIndex:) = 0
          case (ST_INT8)
            rstorageNode%p_Iint8_2D(:,istartIndex:) = 0_I8
          case (ST_INT16)
            rstorageNode%p_Iint16_2D(:,istartIndex:) = 0_I16
          case (ST_INT32)
            rstorageNode%p_Iint32_2D(:,istartIndex:) = 0_I32
          case (ST_INT64)
            rstorageNode%p_Iint64_2D(:,istartIndex:) = 0_I64
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

    case (3)

      select case (cinitNewBlock)
      case (ST_NEWBLOCK_ZERO)
        ! Clear the vector
        if (present(istopIndex)) then
          select case (rstorageNode%idataType)
          case (ST_SINGLE)
            rstorageNode%p_Fsingle3D(:,:,istartIndex:istopIndex) = 0.0_SP
          case (ST_DOUBLE)
            rstorageNode%p_Ddouble3D(:,:,istartIndex:istopIndex) = 0.0_DP
          case (ST_QUAD)
            rstorageNode%p_Qquad3D(:,:,istartIndex:istopIndex) = 0.0_QP
          case (ST_INT)
            rstorageNode%p_Iinteger3D(:,:,istartIndex:istopIndex) = 0
          case (ST_INT8)
            rstorageNode%p_Iint8_3D(:,:,istartIndex:istopIndex) = 0_I8
          case (ST_INT16)
            rstorageNode%p_Iint16_3D(:,:,istartIndex:istopIndex) = 0_I16
          case (ST_INT32)
            rstorageNode%p_Iint32_3D(:,:,istartIndex:istopIndex) = 0_I32
          case (ST_INT64)
            rstorageNode%p_Iint64_3D(:,:,istartIndex:istopIndex) = 0_I64
          case (ST_LOGICAL)
            rstorageNode%p_Blogical3D(:,:,istartIndex:istopIndex) = .false.
          case (ST_CHAR)
            rstorageNode%p_Schar3D(:,:,istartIndex:istopIndex) = achar(0)
          end select
        else
          select case (rstorageNode%idataType)
          case (ST_SINGLE)
            rstorageNode%p_Fsingle3D(:,:,istartIndex:) = 0.0_SP
          case (ST_DOUBLE)
            rstorageNode%p_Ddouble3D(:,:,istartIndex:) = 0.0_DP
          case (ST_QUAD)
            rstorageNode%p_Qquad3D(:,:,istartIndex:) = 0.0_QP
          case (ST_INT)
            rstorageNode%p_Iinteger3D(:,:,istartIndex:) = 0
          case (ST_INT8)
            rstorageNode%p_Iint8_3D(:,:,istartIndex:) = 0_I8
          case (ST_INT16)
            rstorageNode%p_Iint16_3D(:,:,istartIndex:) = 0_I16
          case (ST_INT32)
            rstorageNode%p_Iint32_3D(:,:,istartIndex:) = 0_I32
          case (ST_INT64)
            rstorageNode%p_Iint64_3D(:,:,istartIndex:) = 0_I64
          case (ST_LOGICAL)
            rstorageNode%p_Blogical3D(:,:,istartIndex:) = .false.
          case (ST_CHAR)
            rstorageNode%p_Schar3D(:,:,istartIndex:) = achar(0)
          end select
        end if

      case (ST_NEWBLOCK_ORDERED)
        call output_line ('Ordering not available for multidimensional array!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_initialiseNode')
        call sys_halt()

      end select

    case default
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
  integer, intent(in) :: ihandle

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO, ST_NEWBLOCK_NOINIT,
  ! ST_NEWBLOCK_ORDERED). Specifies how to initialise the data block associated
  ! to ihandle.
  integer, intent(in) :: cinitNewBlock

  ! OPTIONAL: Start index of Block
  integer, intent(in), optional :: istartIndex
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
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
      call storage_initialiseNode (p_rnode,cinitNewBlock,1)
    end if

  end subroutine storage_initialiseBlock

!************************************************************************

!<subroutine>

  subroutine storage_newDefault (scall, sname, Isize, ctype, ihandle, &
                                 cinitNewBlock, rheap)

!<description>
  !This routine reserves a  memory block of desired size and type.
  !It is actually a wrapper for subroutines 'storage_newxD'.
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(in) :: scall

  !clear name of data field
  character(LEN=*), intent(in) :: sname

  !requested storage size
  integer, dimension(:), intent(in) :: Isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(in) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  integer, intent(in) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

!</subroutine>

    select case (size(Isize))
    case (1)
      call storage_new1D (scall, sname, Isize(1), ctype, ihandle, &
                          cinitNewBlock, rheap)
    case (2)
      call storage_new2D (scall, sname, Isize, ctype, ihandle, &
                          cinitNewBlock, rheap)
    case (3)
      call storage_new3D (scall, sname, Isize, ctype, ihandle, &
                          cinitNewBlock, rheap)
    case default
      call output_line ('Memory blocks of dimension larger than 3 is not supported', &
                        OU_CLASS_MSG,OU_MODE_STD,'storage_newDefault')
      ihandle = 0
    end select

  end subroutine storage_newDefault

!************************************************************************

!<subroutine>

  subroutine storage_newDefaultFixed (scall, sname, Ilbound, Iubound, ctype,&
                                      ihandle, cinitNewBlock, rheap)

!<description>
  !This routine reserves a memory block of desired bounds and type.
  !It is actually a wrapper for subroutines 'storage_newxDfixed'.
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(in) :: scall

  !clear name of data field
  character(LEN=*), intent(in) :: sname

  !requested lower bounds for each dimension
  integer, dimension(:), intent(in) :: Ilbound

  !requested upper bounds for each dimension
  integer, dimension(:), intent(in) :: Iubound

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(in) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  integer, intent(in) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

!</subroutine>

    if (size(Ilbound) .ne. size(Iubound)) then
      call output_line ('dim(Ilbound) /= dim(Iubound)!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'storage_newDefaultFixed')
      ihandle = 0
      return
    end if
    
    select case (size(Ilbound))
    case (1)
      call storage_new1Dfixed (scall, sname, Ilbound(1), Iubound(1),&
                               ctype, ihandle, cinitNewBlock, rheap)
    case (2)
      call storage_new2Dfixed (scall, sname, Ilbound, Iubound,&
                               ctype, ihandle, cinitNewBlock, rheap)
    case (3)
      call storage_new3Dfixed (scall, sname, Ilbound, Iubound,&
                               ctype, ihandle, cinitNewBlock, rheap)
    case default
      call output_line ('Memory blocks of dimensions larger than 3 are not supported', &
                        OU_CLASS_WARNING,OU_MODE_STD,'storage_newDefaultFixed')
      ihandle = 0
    end select

  end subroutine storage_newDefaultFixed

!************************************************************************

!<subroutine>

  subroutine storage_new1D (scall, sname, isize, ctype, ihandle, &
                            cinitNewBlock, rheap)

!<description>
  !This routine reserves a 1D memory block of desired size and type.
!</description>

!<input>

  ! name of the calling routine
  character(LEN=*), intent(in) :: scall

  ! clear name of data field
  character(LEN=*), intent(in) :: sname

  ! requested storage size
  integer, intent(in) :: isize

  ! data type, one of the ST_XXXX constants
  integer, intent(in) :: ctype

  ! init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT,ST_NEWBLOCK_ORDERED)
  integer, intent(in) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode
  character(LEN=SYS_NAMELEN) :: snameBackup
  integer :: ier

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
    
#ifdef USE_C_PTR_STORAGE

    ! Allocate memory according to isize:
    select case (ctype)
    case (ST_SINGLE)
      p_rnode%imemBytes = int(isize,I64)*ST_SINGLE_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Fsingle1D, (/isize/))
    case (ST_DOUBLE)
      p_rnode%imemBytes = int(isize,I64)*ST_DOUBLE_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Ddouble1D, (/isize/))
    case (ST_QUAD)
      p_rnode%imemBytes = int(isize,I64)*ST_QUAD_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Qquad1D, (/isize/))
    case (ST_INT)
      p_rnode%imemBytes = int(isize,I64)*ST_INT_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iinteger1D, (/isize/))
    case (ST_INT8)
      p_rnode%imemBytes = int(isize,I64)*ST_INT8_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint8_1D, (/isize/))
    case (ST_INT16)
      p_rnode%imemBytes = int(isize,I64)*ST_INT16_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint16_1D, (/isize/))
    case (ST_INT32)
      p_rnode%imemBytes = int(isize,I64)*ST_INT32_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint32_1D, (/isize/))
    case (ST_INT64)
      p_rnode%imemBytes = int(isize,I64)*ST_INT64_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint64_1D, (/isize/))
    case (ST_LOGICAL)
      p_rnode%imemBytes = int(isize,I64)*ST_LOGICAL_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Blogical1D, (/isize/))
    case (ST_CHAR)
      p_rnode%imemBytes = int(isize,I64)*ST_CHAR_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Schar1D, (/isize/))
    case default
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new1D')
      call sys_halt()
    end select

#else
   
    ! Allocate memory according to isize:
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_SINGLE_BYTES
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_DOUBLE_BYTES
    case (ST_QUAD)
      allocate(p_rnode%p_Qquad1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_QUAD_BYTES
    case (ST_INT)
      allocate(p_rnode%p_Iinteger1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_INT_BYTES
    case (ST_INT8)
      allocate(p_rnode%p_Iint8_1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_INT8_BYTES
    case (ST_INT16)
      allocate(p_rnode%p_Iint16_1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_INT16_BYTES
    case (ST_INT32)
      allocate(p_rnode%p_Iint32_1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_INT32_BYTES
    case (ST_INT64)
      allocate(p_rnode%p_Iint64_1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_INT64_BYTES
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_LOGICAL_BYTES
    case (ST_CHAR)
      allocate(p_rnode%p_Schar1D(isize))
      p_rnode%imemBytes = int(isize,I64)*ST_CHAR_BYTES
    case default
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new1D')
      call sys_halt()
    end select

#endif

    !$omp critical(storage_global_heap_modify)
    p_rheap%itotalMem = p_rheap%itotalMem + p_rnode%imemBytes
    if (p_rheap%itotalMem .gt. p_rheap%itotalMemMax) &
        p_rheap%itotalMemMax = p_rheap%itotalMem
    !$omp end critical(storage_global_heap_modify)
    
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
  character(LEN=*), intent(in) :: scall

  !clear name of data field
  character(LEN=*), intent(in) :: sname

  !requested lower bound
  integer, intent(in) :: ilbound

  !requested upper bound
  integer, intent(in) :: iubound

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(in) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT,ST_NEWBLOCK_ORDERED)
  integer, intent(in) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode
  integer :: isize
  character(LEN=SYS_NAMELEN) :: snameBackup

    ! Can we use the standard routine?
    if (ilbound .eq. 1) then
      call storage_new1D(scall, sname, iubound, ctype, ihandle,&
                         cinitNewBlock, rheap)
      return
    end if
    
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
    
#ifdef USE_C_PTR_STORAGE
    call output_line ('Resorting to standard ALLOCATE for fixed-size memory!', &
                      OU_CLASS_MSG,OU_MODE_STD,'storage_new1Dfixed')
#endif

    ! Allocate memory according to isize:
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_SINGLE_BYTES
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_DOUBLE_BYTES
    case (ST_QUAD)
      allocate(p_rnode%p_Qquad1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_QUAD_BYTES
    case (ST_INT)
      allocate(p_rnode%p_Iinteger1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_INT_BYTES
    case (ST_INT8)
      allocate(p_rnode%p_Iint8_1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_INT8_BYTES
    case (ST_INT16)
      allocate(p_rnode%p_Iint16_1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_INT16_BYTES
    case (ST_INT32)
      allocate(p_rnode%p_Iint32_1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_INT32_BYTES
    case (ST_INT64)
      allocate(p_rnode%p_Iint64_1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_INT64_BYTES
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_LOGICAL_BYTES
    case (ST_CHAR)
      allocate(p_rnode%p_Schar1D(ilbound:iubound))
      p_rnode%imemBytes = int(isize,I64)*ST_CHAR_BYTES
    case default
      call output_line ('Unsupported memory type!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'storage_new1Dfixed')
      call sys_halt()
    end select
    
    !$omp critical(storage_global_heap_modify)
    p_rheap%itotalMem = p_rheap%itotalMem + p_rnode%imemBytes
    if (p_rheap%itotalMem .gt. p_rheap%itotalMemMax) &
        p_rheap%itotalMemMax = p_rheap%itotalMem
    !$omp end critical(storage_global_heap_modify)
    
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
  character(LEN=*), intent(in) :: scall

  !clear name of data field
  character(LEN=*), intent(in) :: sname

  !requested storage size for 1st and 2nd dimension
  integer, dimension(2), intent(in) :: Isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(in) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  integer, intent(in) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode
  character(LEN=SYS_NAMELEN) :: snameBackup
  integer :: ier
    
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
    
#ifdef USE_C_PTR_STORAGE

    ! Allocate memory according to Isize:
    select case (ctype)
    case (ST_SINGLE)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_SINGLE_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Fsingle2D, Isize)
    case (ST_DOUBLE)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_DOUBLE_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Ddouble2D, Isize)
    case (ST_QUAD)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_QUAD_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Qquad2D, Isize)
    case (ST_INT)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iinteger2D, Isize)
    case (ST_INT8)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT8_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint8_2D, Isize)
    case (ST_INT16)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT16_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint16_2D, Isize)
    case (ST_INT32)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT32_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint32_2D, Isize)
    case (ST_INT64)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT64_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint64_2D, Isize)
    case (ST_LOGICAL)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_LOGICAL_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Blogical2D, Isize)
    case (ST_CHAR)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_CHAR_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Schar2D, Isize)
    case default
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new2D')
      call sys_halt()
    end select

#else

    ! Allocate memory according to Isize:
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_SINGLE_BYTES
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_DOUBLE_BYTES
    case (ST_QUAD)
      allocate(p_rnode%p_Qquad2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_QUAD_BYTES
    case (ST_INT)
      allocate(p_rnode%p_Iinteger2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT_BYTES
    case (ST_INT8)
      allocate(p_rnode%p_Iint8_2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT8_BYTES
    case (ST_INT16)
      allocate(p_rnode%p_Iint16_2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT16_BYTES
    case (ST_INT32)
      allocate(p_rnode%p_Iint32_2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT32_BYTES
    case (ST_INT64)
      allocate(p_rnode%p_Iint64_2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT64_BYTES
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_LOGICAL_BYTES
    case (ST_CHAR)
      allocate(p_rnode%p_Schar2D(Isize(1),Isize(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_CHAR_BYTES
    case default
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new2D')
      call sys_halt()
    end select

#endif
    
    !$omp critical(storage_global_heap_modify)
    p_rheap%itotalMem = p_rheap%itotalMem + p_rnode%imemBytes
    if (p_rheap%itotalMem .gt. p_rheap%itotalMemMax) &
        p_rheap%itotalMemMax = p_rheap%itotalMem
    !$omp end critical(storage_global_heap_modify)
    
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
  character(LEN=*), intent(in) :: scall

  !clear name of data field
  character(LEN=*), intent(in) :: sname

  !requested lower bounds for 1st and 2nd dimension
  integer, dimension(2), intent(in) :: Ilbound

  !requested upper bounds for 1st and 2nd dimension
  integer, dimension(2), intent(in) :: Iubound

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(in) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  integer, intent(in) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

!</subroutine>

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    integer, dimension(2) :: Isize
    character(LEN=SYS_NAMELEN) :: snameBackup
    
    ! Can we use the standard routine?
    if (all(Ilbound .eq. 1)) then
      call storage_new2D(scall, sname, Iubound, ctype, ihandle,&
                         cinitNewBlock, rheap)
      return
    end if

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
        
#ifdef USE_C_PTR_STORAGE
    call output_line ('Resorting to standard ALLOCATE for fixed-size memory!', &
                      OU_CLASS_MSG,OU_MODE_STD,'storage_new2Dfixed')
#endif

    ! Allocate memory according to Isize:
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_SINGLE_BYTES
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_DOUBLE_BYTES
    case (ST_QUAD)
      allocate(p_rnode%p_Qquad2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_QUAD_BYTES
    case (ST_INT)
      allocate(p_rnode%p_Iinteger2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT_BYTES
    case (ST_INT8)
      allocate(p_rnode%p_Iint8_2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT8_BYTES
    case (ST_INT16)
      allocate(p_rnode%p_Iint16_2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT16_BYTES
    case (ST_INT32)
      allocate(p_rnode%p_Iint32_2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT32_BYTES
    case (ST_INT64)
      allocate(p_rnode%p_Iint64_2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_INT64_BYTES
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_LOGICAL_BYTES
    case (ST_CHAR)
      allocate(p_rnode%p_Schar2D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*ST_CHAR_BYTES
    case default
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new2Dfixed')
      call sys_halt()
    end select
    
    !$omp critical(storage_global_heap_modify)
    p_rheap%itotalMem = p_rheap%itotalMem + p_rnode%imemBytes
    if (p_rheap%itotalMem .gt. p_rheap%itotalMemMax) &
        p_rheap%itotalMemMax = p_rheap%itotalMem
    !$omp end critical(storage_global_heap_modify)
    
    ! Initialise the storage block
    call storage_initialiseBlock (ihandle, cinitNewBlock, rheap, Ilbound(2))
    
  end subroutine storage_new2Dfixed

!************************************************************************

!<subroutine>

  subroutine storage_new3D (scall, sname, Isize, ctype, ihandle, &
                            cinitNewBlock, rheap)

!<description>
  !This routine reserves a 3D memory block of desired size and type.
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(in) :: scall

  !clear name of data field
  character(LEN=*), intent(in) :: sname

  !requested storage size for 1st, 2nd and 3rd dimension
  integer, dimension(3), intent(in) :: Isize

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(in) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  integer, intent(in) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode
  character(LEN=SYS_NAMELEN) :: snameBackup
  integer :: ier

    if ((Isize(1) .eq. 0) .or. (Isize(2) .eq. 0) .or. (Isize(3) .eq. 0)) then
      call output_line ('Isize=0!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'storage_new3D')
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
    p_rnode%idimension = 3
    p_rnode%sname = snameBackup
    
#ifdef USE_C_PTR_STORAGE

    ! Allocate memory according to Isize:
    select case (ctype)
    case (ST_SINGLE)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_SINGLE_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Fsingle3D, Isize)
    case (ST_DOUBLE)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_DOUBLE_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Ddouble3D, Isize)
    case (ST_QUAD)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_QUAD_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Qquad3D, Isize)
    case (ST_INT)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iinteger3D, Isize)
    case (ST_INT8)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT8_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint8_3D, Isize)
    case (ST_INT16)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT16_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint16_3D, Isize)
    case (ST_INT32)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT32_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint32_3D, Isize)
    case (ST_INT64)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT64_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Iint64_3D, Isize)
    case (ST_LOGICAL)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_LOGICAL_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Blogical3D, Isize)
    case (ST_CHAR)
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_CHAR_BYTES
      ier = c_allocate(p_rnode%chostMemPtr, p_rnode%imemBytes)
      call c_f_pointer(p_rnode%chostMemPtr, p_rnode%p_Schar3D, Isize)
    case default
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new3D')
      call sys_halt()
    end select

#else

    ! Allocate memory according to Isize:
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_SINGLE_BYTES
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_DOUBLE_BYTES
    case (ST_QUAD)
      allocate(p_rnode%p_Qquad3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_QUAD_BYTES
    case (ST_INT)
      allocate(p_rnode%p_Iinteger3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT_BYTES
    case (ST_INT8)
      allocate(p_rnode%p_Iint8_3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT8_BYTES
    case (ST_INT16)
      allocate(p_rnode%p_Iint16_3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT16_BYTES
    case (ST_INT32)
      allocate(p_rnode%p_Iint32_3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT32_BYTES
    case (ST_INT64)
      allocate(p_rnode%p_Iint64_3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT64_BYTES
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_LOGICAL_BYTES
    case (ST_CHAR)
      allocate(p_rnode%p_Schar3D(Isize(1),Isize(2),Isize(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_CHAR_BYTES
    case default
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new3D')
      call sys_halt()
    end select

#endif
    
    !$omp critical(storage_global_heap_modify)
    p_rheap%itotalMem = p_rheap%itotalMem + p_rnode%imemBytes
    if (p_rheap%itotalMem .gt. p_rheap%itotalMemMax) &
        p_rheap%itotalMemMax = p_rheap%itotalMem
    !$omp end critical(storage_global_heap_modify)
    
    ! Initialise the storage block
    call storage_initialiseBlock (ihandle, cinitNewBlock, rheap)
    
  end subroutine storage_new3D

!************************************************************************

!<subroutine>

  subroutine storage_new3Dfixed (scall, sname, Ilbound, Iubound, ctype,&
                                 ihandle, cinitNewBlock, rheap)

!<description>
  !This routine reserves a 3D memory block of desired bounds and type.
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(in) :: scall

  !clear name of data field
  character(LEN=*), intent(in) :: sname

  !requested lower bounds for 1st, 2nd and 3rd dimension
  integer, dimension(3), intent(in) :: Ilbound

  !requested upper bounds for 1st, 2nd and 3rd dimension
  integer, dimension(3), intent(in) :: Iubound

  !data type (ST_SINGLE,ST_DOUBLE,ST_INT,ST_LOGICAL,ST_CHAR)
  integer, intent(in) :: ctype

  !init new storage block (ST_NEWBLOCK_ZERO,ST_NEWBLOCK_NOINIT)
  integer, intent(in) :: cinitNewBlock

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

!</subroutine>

    ! Pointer to the heap
    type(t_storageBlock), pointer :: p_rheap
    type(t_storageNode), pointer :: p_rnode
    integer, dimension(3) :: Isize
    character(LEN=SYS_NAMELEN) :: snameBackup
    
    ! Can we use the standard routine?
    if (all(Ilbound .eq. 1)) then
      call storage_new3D(scall, sname, Iubound, ctype, ihandle,&
                         cinitNewBlock, rheap)
      return
    end if

    Isize=Iubound-Ilbound+1
    if ((Isize(1) .eq. 0) .or. (Isize(2) .eq. 0) .or. (Isize(3) .eq. 0)) then
      call output_line ('Isize=0!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'storage_new3Dfixed')
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
    p_rnode%idimension = 3
    p_rnode%sname = snameBackup
    
#ifdef USE_C_PTR_STORAGE
    call output_line ('Resorting to standard ALLOCATE for fixed-size memory!', &
                      OU_CLASS_MSG,OU_MODE_STD,'storage_new3Dfixed')
#endif

    ! Allocate memory according to Isize:
    select case (ctype)
    case (ST_SINGLE)
      allocate(p_rnode%p_Fsingle3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                   Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_SINGLE_BYTES
    case (ST_DOUBLE)
      allocate(p_rnode%p_Ddouble3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                   Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_DOUBLE_BYTES
    case (ST_QUAD)
      allocate(p_rnode%p_Qquad3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                 Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_QUAD_BYTES
    case (ST_INT)
      allocate(p_rnode%p_Iinteger3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                    Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT_BYTES
    case (ST_INT8)
      allocate(p_rnode%p_Iint8_3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                  Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT8_BYTES
    case (ST_INT16)
      allocate(p_rnode%p_Iint16_3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                   Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT16_BYTES
    case (ST_INT32)
      allocate(p_rnode%p_Iint32_3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                   Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT32_BYTES
    case (ST_INT64)
      allocate(p_rnode%p_Iint64_3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                   Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_INT64_BYTES
    case (ST_LOGICAL)
      allocate(p_rnode%p_Blogical3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                    Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_LOGICAL_BYTES
    case (ST_CHAR)
      allocate(p_rnode%p_Schar3D(Ilbound(1):Iubound(1),Ilbound(2):Iubound(2),&
                                 Ilbound(3):Iubound(3)))
      p_rnode%imemBytes = int(Isize(1),I64)*int(Isize(2),I64)*&
                          int(Isize(3),I64)*ST_CHAR_BYTES
    case default
      call output_line ('Unsupported memory type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_new2Dfixed')
      call sys_halt()
    end select
    
    !$omp critical(storage_global_heap_modify)
    p_rheap%itotalMem = p_rheap%itotalMem + p_rnode%imemBytes
    if (p_rheap%itotalMem .gt. p_rheap%itotalMemMax) &
        p_rheap%itotalMemMax = p_rheap%itotalMem
    !$omp end critical(storage_global_heap_modify)
    
    ! Initialise the storage block
    call storage_initialiseBlock (ihandle, cinitNewBlock, rheap, Ilbound(2))
    
  end subroutine storage_new3Dfixed

!************************************************************************

!<subroutine>

  subroutine storage_newIndirect (scall, sname, ihandleTemplate, ihandle, rheap)

!<description>
  !This routine reserves a memory block of desired size and type
  !adopting all dimensions, types, etc. from the template handle
!</description>

!<input>

  !name of the calling routine
  character(LEN=*), intent(in) :: scall

  !clear name of data field
  character(LEN=*), intent(in) :: sname

  ! Handle of the template memory block.
  integer, intent(in) :: ihandleTemplate

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</inputoutput>

!<output>

  ! Handle of the memory block.
  integer, intent(out) :: ihandle

!</output>

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

    if (ihandleTemplate .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_newIndirect')
      call sys_halt()
    end if
    if (.not. associated(p_rheap%p_Rdescriptors)) then
      call output_line ('Heap not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_newIndirect')
      call sys_halt()
    end if

    p_rnode => p_rheap%p_Rdescriptors(ihandleTemplate)

    ! What dimension are we?
    select case(p_rnode%idimension)

    case (1)
      
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Fsingle1D,1),&
                          ubound(p_rnode%p_Fsingle1D,1),&
                          ST_SINGLE, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_DOUBLE)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Ddouble1D,1),&
                          ubound(p_rnode%p_Ddouble1D,1),&
                          ST_DOUBLE, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_QUAD)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Qquad1D,1),&
                          ubound(p_rnode%p_Qquad1D,1),&
                          ST_QUAD, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Iinteger1D,1),&
                          ubound(p_rnode%p_Iinteger1D,1),&
                          ST_INT, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT8)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Iint8_1D,1),&
                          ubound(p_rnode%p_Iint8_1D,1),&
                          ST_INT8, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT16)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Iint16_1D,1),&
                          ubound(p_rnode%p_Iint16_1D,1),&
                          ST_INT16, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT32)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Iint32_1D,1),&
                          ubound(p_rnode%p_Iint32_1D,1),&
                          ST_INT32, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT64)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Iint64_1D,1),&
                          ubound(p_rnode%p_Iint64_1D,1),&
                          ST_INT64, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_LOGICAL)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Blogical1D,1),&
                          ubound(p_rnode%p_Blogical1D,1),&
                          ST_LOGICAL, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_CHAR)
        call storage_new ('storage_newIndirect',sname,&
                          lbound(p_rnode%p_Schar1D,1),&
                          ubound(p_rnode%p_Schar1D,1),&
                          ST_CHAR, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      end select

    case (2)
      
      select case (p_rnode%IdataType)
      case (ST_SINGLE)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Fsingle2D),&
                          ubound(p_rnode%p_Fsingle2D),&
                          ST_SINGLE, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_DOUBLE)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Ddouble2D),&
                          ubound(p_rnode%p_Ddouble2D),&
                          ST_DOUBLE, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_QUAD)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Qquad2D),&
                          ubound(p_rnode%p_Qquad2D),&
                          ST_QUAD, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iinteger2D),&
                          ubound(p_rnode%p_Iinteger2D),&
                          ST_INT, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT8)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iint8_2D),&
                          ubound(p_rnode%p_Iint8_2D),&
                          ST_INT8, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT16)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iint16_2D),&
                          ubound(p_rnode%p_Iint16_2D),&
                          ST_INT16, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT32)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iint32_2D),&
                          ubound(p_rnode%p_Iint32_2D),&
                          ST_INT32, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT64)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iint64_2D),&
                          ubound(p_rnode%p_Iint64_2D),&
                          ST_INT64, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_LOGICAL)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Blogical2D),&
                          ubound(p_rnode%p_Blogical2D),&
                          ST_LOGICAL, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_CHAR)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Schar2D),&
                          ubound(p_rnode%p_Schar2D),&
                          ST_CHAR, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      end select

    case (3)

      select case (p_rnode%IdataType)
      case (ST_SINGLE)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Fsingle3D),&
                          ubound(p_rnode%p_Fsingle3D),&
                          ST_SINGLE, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_DOUBLE)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Ddouble3D),&
                          ubound(p_rnode%p_Ddouble3D),&
                          ST_DOUBLE, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_QUAD)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Qquad3D),&
                          ubound(p_rnode%p_Qquad3D),&
                          ST_QUAD, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iinteger3D),&
                          ubound(p_rnode%p_Iinteger3D),&
                          ST_INT, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT8)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iint8_3D),&
                          ubound(p_rnode%p_Iint8_3D),&
                          ST_INT8, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT16)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iint16_3D),&
                          ubound(p_rnode%p_Iint16_3D),&
                          ST_INT16, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT32)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iint32_3D),&
                          ubound(p_rnode%p_Iint32_3D),&
                          ST_INT32, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT64)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Iint64_3D),&
                          ubound(p_rnode%p_Iint64_3D),&
                          ST_INT64, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_LOGICAL)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Blogical3D),&
                          ubound(p_rnode%p_Blogical3D),&
                          ST_LOGICAL, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_CHAR)
        call storage_new ('storage_newIndirect', sname,&
                          lbound(p_rnode%p_Schar3D),&
                          ubound(p_rnode%p_Schar3D),&
                          ST_CHAR, ihandle, ST_NEWBLOCK_NOINIT, p_rheap)
      end select

    case default
      call output_line ('Unsupported dimension!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'storage_newIndirect')
      call sys_halt()
    end select

  end subroutine storage_newIndirect

!************************************************************************

!<subroutine>

  subroutine storage_free (ihandle, rheap)

!<description>
  ! This routine releases a handle from a heap and deallocates the
  ! associated memory. ihandle is set to ST_NOHANDLE upon return.
!</description>

!<inputoutput>

  ! Handle of the memory block to be releases
  integer, intent(inout) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

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
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Releasing ST_NOHANDLE is not allowed!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_free')
      call sys_halt()
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Is the node associated at all?
    if (p_rnode%idataType .le. ST_NOHANDLE) then
      call output_line ('Trying to release nonexistent handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_free')
      call output_line ('Handle number: '//trim(sys_siL(ihandle,11)), &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_free')
      call sys_halt()
    end if

#ifdef USE_C_PTR_STORAGE
    ! Release host memory physically
    if (storage_isAssociated(p_rnode%chostMemPtr)) then
      if (c_deallocate(p_rnode%chostMemPtr) .gt. 0) then
        call output_line ('Error in freeing memory with c_deallocate!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_free')
        call sys_halt()
      end if
    end if

    ! Nullify pointers
    nullify(p_rnode%p_Fsingle1D)
    nullify(p_rnode%p_Ddouble1D)
    nullify(p_rnode%p_Qquad1D)
    nullify(p_rnode%p_Iinteger1D)
    nullify(p_rnode%p_Iint8_1D)
    nullify(p_rnode%p_Iint16_1D)
    nullify(p_rnode%p_Iint32_1D)
    nullify(p_rnode%p_Iint64_1D)
    nullify(p_rnode%p_Blogical1D)
    nullify(p_rnode%p_Schar1D)

    nullify(p_rnode%p_Fsingle2D)
    nullify(p_rnode%p_Ddouble2D)
    nullify(p_rnode%p_Qquad2D)
    nullify(p_rnode%p_Iinteger2D)
    nullify(p_rnode%p_Iint8_2D)
    nullify(p_rnode%p_Iint16_2D)
    nullify(p_rnode%p_Iint32_2D)
    nullify(p_rnode%p_Iint64_2D)
    nullify(p_rnode%p_Blogical2D)
    nullify(p_rnode%p_Schar2D)

    nullify(p_rnode%p_Fsingle3D)
    nullify(p_rnode%p_Ddouble3D)
    nullify(p_rnode%p_Qquad3D)
    nullify(p_rnode%p_Iinteger3D)
    nullify(p_rnode%p_Iint8_3D)
    nullify(p_rnode%p_Iint16_3D)
    nullify(p_rnode%p_Iint32_3D)
    nullify(p_rnode%p_Iint64_3D)
    nullify(p_rnode%p_Blogical3D)
    nullify(p_rnode%p_Schar3D)
#endif

#ifdef ENABLE_COPROCESSOR_SUPPORT
    if (storage_isAssociated(p_rnode%cdeviceMemPtr))&
        call coproc_freeMemoryOnDevice(p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
#endif

    ! Release the memory assigned to that
    if (associated(p_rnode%p_Fsingle1D))  deallocate(p_rnode%p_Fsingle1D)
    if (associated(p_rnode%p_Ddouble1D))  deallocate(p_rnode%p_Ddouble1D)
    if (associated(p_rnode%p_Qquad1D))    deallocate(p_rnode%p_Qquad1D)
    if (associated(p_rnode%p_Iinteger1D)) deallocate(p_rnode%p_Iinteger1D)
    if (associated(p_rnode%p_Iint8_1D))   deallocate(p_rnode%p_Iint8_1D)
    if (associated(p_rnode%p_Iint16_1D))  deallocate(p_rnode%p_Iint16_1D)
    if (associated(p_rnode%p_Iint32_1D))  deallocate(p_rnode%p_Iint32_1D)
    if (associated(p_rnode%p_Iint64_1D))  deallocate(p_rnode%p_Iint64_1D)
    if (associated(p_rnode%p_Blogical1D)) deallocate(p_rnode%p_Blogical1D)
    if (associated(p_rnode%p_Schar1D))    deallocate(p_rnode%p_Schar1D)

    if (associated(p_rnode%p_Fsingle2D))  deallocate(p_rnode%p_Fsingle2D)
    if (associated(p_rnode%p_Ddouble2D))  deallocate(p_rnode%p_Ddouble2D)
    if (associated(p_rnode%p_Qquad2D))    deallocate(p_rnode%p_Qquad2D)
    if (associated(p_rnode%p_Iinteger2D)) deallocate(p_rnode%p_Iinteger2D)
    if (associated(p_rnode%p_Iint8_2D))   deallocate(p_rnode%p_Iint8_2D)
    if (associated(p_rnode%p_Iint16_2D))  deallocate(p_rnode%p_Iint16_2D)
    if (associated(p_rnode%p_Iint32_2D))  deallocate(p_rnode%p_Iint32_2D)
    if (associated(p_rnode%p_Iint64_2D))  deallocate(p_rnode%p_Iint64_2D)
    if (associated(p_rnode%p_Blogical2D)) deallocate(p_rnode%p_Blogical2D)
    if (associated(p_rnode%p_Schar2D))    deallocate(p_rnode%p_Schar2D)

    if (associated(p_rnode%p_Fsingle3D))  deallocate(p_rnode%p_Fsingle3D)
    if (associated(p_rnode%p_Ddouble3D))  deallocate(p_rnode%p_Ddouble3D)
    if (associated(p_rnode%p_Qquad3D))    deallocate(p_rnode%p_Qquad3D)
    if (associated(p_rnode%p_Iinteger3D)) deallocate(p_rnode%p_Iinteger3D)
    if (associated(p_rnode%p_Iint8_3D))   deallocate(p_rnode%p_Iint8_3D)
    if (associated(p_rnode%p_Iint16_3D))  deallocate(p_rnode%p_Iint16_3D)
    if (associated(p_rnode%p_Iint32_3D))  deallocate(p_rnode%p_Iint32_3D)
    if (associated(p_rnode%p_Iint64_3D))  deallocate(p_rnode%p_Iint64_3D)
    if (associated(p_rnode%p_Blogical3D)) deallocate(p_rnode%p_Blogical3D)
    if (associated(p_rnode%p_Schar3D))    deallocate(p_rnode%p_Schar3D)
    
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

!<input>

    ! Handle of the memory block to be cleared
  integer, intent(in) :: ihandle

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

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
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'storage_clean')
      call sys_halt()
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Is the node associated at all?
    if (p_rnode%idataType .le. ST_NOHANDLE) then
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
      case (ST_QUAD)
        p_rnode%p_Qquad1D = 0.0_QP
      case (ST_INT)
        p_rnode%p_Iinteger1D = 0
      case (ST_INT8)
        p_rnode%p_Iint8_1D = 0_I8
      case (ST_INT16)
        p_rnode%p_Iint16_1D = 0_I16
      case (ST_INT32)
        p_rnode%p_Iint32_1D = 0_I32
      case (ST_INT64)
        p_rnode%p_Iint64_1D = 0_I64
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
      case (ST_QUAD)
        p_rnode%p_Qquad2D = 0.0_QP
      case (ST_INT)
        p_rnode%p_Iinteger2D = 0
      case (ST_INT8)
        p_rnode%p_Iint8_2D = 0_I8
      case (ST_INT16)
        p_rnode%p_Iint16_2D = 0_I16
      case (ST_INT32)
        p_rnode%p_Iint32_2D = 0_I32
      case (ST_INT64)
        p_rnode%p_Iint64_2D = 0_I64
      case (ST_LOGICAL)
        p_rnode%p_Blogical2D = .false.
      case (ST_CHAR)
        p_rnode%p_Schar2D = achar(0)
      end select

    case (3)
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        p_rnode%p_Fsingle3D = 0.0_SP
      case (ST_DOUBLE)
        p_rnode%p_Ddouble3D = 0.0_DP
      case (ST_QUAD)
        p_rnode%p_Qquad3D = 0.0_QP
      case (ST_INT)
        p_rnode%p_Iinteger3D = 0
      case (ST_INT8)
        p_rnode%p_Iint8_3D = 0_I8
      case (ST_INT16)
        p_rnode%p_Iint16_3D = 0_I16
      case (ST_INT32)
        p_rnode%p_Iint32_3D = 0_I32
      case (ST_INT64)
        p_rnode%p_Iint64_3D = 0_I64
      case (ST_LOGICAL)
        p_rnode%p_Blogical3D = .false.
      case (ST_CHAR)
        p_rnode%p_Schar3D = achar(0)
      end select

    case default
      call output_line ('Invalid dimension!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_clear')
      call sys_halt()
    end select
    
  end subroutine storage_clear

!************************************************************************

!<subroutine>

  subroutine storage_getsize1D (ihandle, isize, rheap)

!<description>
  ! Returns the length of a 1-dimensional array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<output>
  ! Length of the array identified by ihandle.
  integer, intent(out) :: isize
!</output>

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
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsize1D')
      call sys_halt()
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Is the node associated at all?
    if (p_rnode%idataType .le. ST_NOHANDLE) then
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
    case (ST_QUAD)
      isize = size(p_rnode%p_Qquad1D)
    case (ST_INT)
      isize = size(p_rnode%p_Iinteger1D)
    case (ST_INT8)
      isize = size(p_rnode%p_Iint8_1D)
    case (ST_INT16)
      isize = size(p_rnode%p_Iint16_1D)
    case (ST_INT32)
      isize = size(p_rnode%p_Iint32_1D)
    case (ST_INT64)
      isize = size(p_rnode%p_Iint64_1D)
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

  subroutine storage_getsizeDefault (ihandle, Isize, rheap)

!<description>
  ! Returns the length of a multidimensional array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<output>
  ! Length of each dimension of the array identified by ihandle.
  ! Note that if the dimension of Isize is not equal to the
  ! dimension of the array an error is thrown.
  integer, dimension(:), intent(out) :: Isize
!</output>

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
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
      call sys_halt()
    end if
    
    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Is the node associated at all?
    if (p_rnode%idataType .le. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
      call output_line ('Handle number: '//trim(sys_siL(ihandle,11)), &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
      call sys_halt()
    end if
    
    ! What are we?
    select case (p_rnode%idimension)
    case (1)
      if (size(Isize) .ne. 1) then
        call output_line ('dim(Isize) /= 1!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
        call sys_halt()
      end if

      select case (p_rnode%idataType)
      case (ST_SINGLE)
        Isize = shape(p_rnode%p_Fsingle1D)
      case (ST_DOUBLE)
        Isize = shape(p_rnode%p_Ddouble1D)
      case (ST_QUAD)
        Isize = shape(p_rnode%p_Qquad1D)
      case (ST_INT)
        Isize = shape(p_rnode%p_Iinteger1D)
      case (ST_INT8)
        Isize = shape(p_rnode%p_Iint8_1D)
      case (ST_INT16)
        Isize = shape(p_rnode%p_Iint16_1D)
      case (ST_INT32)
        Isize = shape(p_rnode%p_Iint32_1D)
      case (ST_INT64)
        Isize = shape(p_rnode%p_Iint64_1D)
      case (ST_LOGICAL)
        Isize = shape(p_rnode%p_Blogical1D)
      case (ST_CHAR)
        Isize = shape(p_rnode%p_Schar1D)
      case default
        call output_line ('Invalid data type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
        call sys_halt()
      end select
      
    case (2)
      if (size(Isize) .ne. 2) then
        call output_line ('dim(Isize) /= 2!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
        call sys_halt()
      end if

      select case (p_rnode%idataType)
      case (ST_SINGLE)
        Isize = shape(p_rnode%p_Fsingle2D)
      case (ST_DOUBLE)
        Isize = shape(p_rnode%p_Ddouble2D)
      case (ST_QUAD)
        Isize = shape(p_rnode%p_Qquad2D)
      case (ST_INT)
        Isize = shape(p_rnode%p_Iinteger2D)
      case (ST_INT8)
        Isize = shape(p_rnode%p_Iint8_2D)
      case (ST_INT16)
        Isize = shape(p_rnode%p_Iint16_2D)
      case (ST_INT32)
        Isize = shape(p_rnode%p_Iint32_2D)
      case (ST_INT64)
        Isize = shape(p_rnode%p_Iint64_2D)
      case (ST_LOGICAL)
        Isize = shape(p_rnode%p_Blogical2D)
      case (ST_CHAR)
        Isize = shape(p_rnode%p_Schar2D)
      case default
        call output_line ('Invalid data type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
        call sys_halt()
      end select
      
    case (3)
      if (size(Isize) .ne. 3) then
        call output_line ('dim(Isize) /= 3!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
        call sys_halt()
      end if

      select case (p_rnode%idataType)
      case (ST_SINGLE)
        Isize = shape(p_rnode%p_Fsingle3D)
      case (ST_DOUBLE)
        Isize = shape(p_rnode%p_Ddouble3D)
      case (ST_QUAD)
        Isize = shape(p_rnode%p_Qquad3D)
      case (ST_INT)
        Isize = shape(p_rnode%p_Iinteger3D)
      case (ST_INT8)
        Isize = shape(p_rnode%p_Iint8_3D)
      case (ST_INT16)
        Isize = shape(p_rnode%p_Iint16_3D)
      case (ST_INT32)
        Isize = shape(p_rnode%p_Iint32_3D)
      case (ST_INT64)
        Isize = shape(p_rnode%p_Iint64_3D)
      case (ST_LOGICAL)
        Isize = shape(p_rnode%p_Blogical3D)
      case (ST_CHAR)
        Isize = shape(p_rnode%p_Schar3D)
      case default
        call output_line ('Invalid data type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
        call sys_halt()
      end select

    case default
      call output_line ('Handle '//trim(sys_siL(ihandle,11))//' is not a valid handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getsizeDefault')
      call sys_halt()
    end select
        
  end subroutine storage_getsizeDefault

!************************************************************************

!<subroutine>

  subroutine storage_getbase_intDefault (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
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

  end subroutine storage_getbase_intDefault

!************************************************************************

!<subroutine>

  subroutine storage_getbase_intUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
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
  ! integer array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle
  
  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
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
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_intUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D(lbnd:ubnd)

  end subroutine storage_getbase_intLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_singleDefault (ihandle, p_Sarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase__singleDefault')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_singleDefault')
      call sys_halt()
    end if

    ! Get the pointer
    p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D

  end subroutine storage_getbase_singleDefault

!************************************************************************

!<subroutine>

  subroutine storage_getbase_singleUBnd (ihandle, p_Sarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
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
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
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
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_singleUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D(lbnd:ubnd)

  end subroutine storage_getbase_singleLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_doubleDefault (ihandle, p_Darray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(DP), dimension(:), pointer :: p_Darray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleDefault')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleDefault')
      call sys_halt()
    end if

    ! Get the pointer
    p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D

  end subroutine storage_getbase_doubleDefault

!************************************************************************

!<subroutine>

  subroutine storage_getbase_doubleUBnd (ihandle, p_Darray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(DP), dimension(:), pointer :: p_Darray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
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
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(DP), dimension(:), pointer :: p_Darray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
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
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_doubleUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D(lbnd:ubnd)

  end subroutine storage_getbase_doubleLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quadDefault (ihandle, p_Qarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(QP), dimension(:), pointer :: p_Qarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quadDefault')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quadDefault')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad1D

  end subroutine storage_getbase_quadDefault

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quadUBnd (ihandle, p_Qarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(QP), dimension(:), pointer :: p_Qarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quadUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quadUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Qquad1D))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quadUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad1D(:ubnd)

  end subroutine storage_getbase_quadUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quadLUBnd (ihandle, p_Qarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(QP), dimension(:), pointer :: p_Qarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quadLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quadLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Qquad1D)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Qquad1D)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quadLUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad1D(lbnd:ubnd)

  end subroutine storage_getbase_quadLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logicalDefault (ihandle, p_Larray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logicalDefault')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logicalDefault')
      call sys_halt()
    end if

    ! Get the pointer
    p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D

  end subroutine storage_getbase_logicalDefault

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logicalUBnd (ihandle, p_Larray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logicalUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logicalUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logicalUBnd')
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
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logicalLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logicalLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logicalUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D(lbnd:ubnd)

  end subroutine storage_getbase_logicalLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_charDefault (ihandle, p_Carray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_charDefault')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_charDefault')
      call sys_halt()
    end if

    ! Get the pointer
    p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D

  end subroutine storage_getbase_charDefault

!************************************************************************

!<subroutine>

  subroutine storage_getbase_charUBnd (ihandle, p_Carray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle
  
  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_charUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_charUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_charUBnd')
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
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_charLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_charLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_charUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D(lbnd:ubnd)

  end subroutine storage_getbase_charLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8Default (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8Default')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8Default')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_1D

  end subroutine storage_getbase_int8Default

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8UBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8UBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8UBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_1D))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8UBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_1D(:ubnd)

  end subroutine storage_getbase_int8UBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8LUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle
  
  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8LUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8LUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_1D)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_1D)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8UBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_1D(lbnd:ubnd)

  end subroutine storage_getbase_int8LUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16Default (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16Default')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16Default')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_1D

  end subroutine storage_getbase_int16Default

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16UBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16UBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16UBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_1D))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16UBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_1D(:ubnd)

  end subroutine storage_getbase_int16UBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16LUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle
  
  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16LUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16LUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_1D)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_1D)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16UBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_1D(lbnd:ubnd)

  end subroutine storage_getbase_int16LUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32Default (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32Default')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32Default')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_1D

  end subroutine storage_getbase_int32Default

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32UBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32UBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32UBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_1D))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32UBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_1D(:ubnd)

  end subroutine storage_getbase_int32UBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32LUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle
  
  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32LUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32LUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_1D)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_1D)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32UBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_1D(lbnd:ubnd)

  end subroutine storage_getbase_int32LUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64Default (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64Default')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64Default')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_1D

  end subroutine storage_getbase_int64Default

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64UBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64UBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64UBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_1D))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64UBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_1D(:ubnd)

  end subroutine storage_getbase_int64UBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64LUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bound.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle
  
  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64LUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64LUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_1D)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. size(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_1D)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64UBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_1D(lbnd:ubnd)

  end subroutine storage_getbase_int64LUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int2DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D

  end subroutine storage_getbase_int2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int2DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int2DUBnd')
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
  ! integer array and adopts the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D(:,lbnd:ubnd)

  end subroutine storage_getbase_int2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single2DDef (ihandle, p_Sarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! single precision array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D

  end subroutine storage_getbase_single2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single2DUBnd (ihandle, p_Sarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single precision array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single2DUBnd')
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
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd
  
  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D(:,lbnd:ubnd)

  end subroutine storage_getbase_single2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double2DDef (ihandle, p_Darray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D

  end subroutine storage_getbase_double2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double2DUBnd (ihandle, p_Darray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array and adopt the given upper bound for the second dimension.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2))) then
       call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double2DUBnd')
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
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D(:,lbnd:ubnd)

  end subroutine storage_getbase_double2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quad2DDef (ihandle, p_Qarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad precision array.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(QP), dimension(:,:), pointer :: p_Qarray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad2D

  end subroutine storage_getbase_quad2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quad2DUBnd (ihandle, p_Qarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad precision array and adopt the given upper bound for the second dimension.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(QP), dimension(:,:), pointer :: p_Qarray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Qquad2D,2))) then
       call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad2D(:,:ubnd)

  end subroutine storage_getbase_quad2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quad2DLUBnd (ihandle, p_Qarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad precision array and adopt the given lower and upper bounds
  ! for the second dimension.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(QP), dimension(:,:), pointer :: p_Qarray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Qquad2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Qquad2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad2D(:,lbnd:ubnd)

  end subroutine storage_getbase_quad2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical2DDef (ihandle, p_Larray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D

  end subroutine storage_getbase_logical2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical2DUBnd (ihandle, p_Larray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D,2))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical2DUBnd')
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
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D(:,lbnd:ubnd)

  end subroutine storage_getbase_logical2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char2DDef (ihandle, p_Carray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D

  end subroutine storage_getbase_char2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char2DUBnd (ihandle, p_Carray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopt the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Schar2D,2))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char2DUBnd')
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
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Schar2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Schar2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D(:,lbnd:ubnd)

  end subroutine storage_getbase_char2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8_2DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_2D

  end subroutine storage_getbase_int8_2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8_2DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_2D,2))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_2D(:,:ubnd)

  end subroutine storage_getbase_int8_2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8_2DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_2D(:,lbnd:ubnd)

  end subroutine storage_getbase_int8_2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16_2DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_2D

  end subroutine storage_getbase_int16_2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16_2DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_2D,2))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_2D(:,:ubnd)

  end subroutine storage_getbase_int16_2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16_2DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_2D(:,lbnd:ubnd)

  end subroutine storage_getbase_int16_2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32_2DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_2D

  end subroutine storage_getbase_int32_2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32_2DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_2D,2))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_2D(:,:ubnd)

  end subroutine storage_getbase_int32_2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32_2DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_2D(:,lbnd:ubnd)

  end subroutine storage_getbase_int32_2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64_2DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_2DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_2DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_2D

  end subroutine storage_getbase_int64_2DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64_2DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_2DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_2DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_2D,2))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_2D(:,:ubnd)

  end subroutine storage_getbase_int64_2DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64_2DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the second dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_2DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_2DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_2D,2)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_2D,2)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_2DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_2D(:,lbnd:ubnd)

  end subroutine storage_getbase_int64_2DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int3DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger3D

  end subroutine storage_getbase_int3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int3DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger3D,3))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger3D(:,:,:ubnd)

  end subroutine storage_getbase_int3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int3DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer, dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iinteger3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_int3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single3DDef (ihandle, p_Sarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! single precision array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle3D

  end subroutine storage_getbase_single3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single3DUBnd (ihandle, p_Sarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! single precision array and adopt the given upper bound for the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle3D,3))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle3D(:,:,:ubnd)

  end subroutine storage_getbase_single3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_single3DLUBnd (ihandle, p_Sarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! single precision array and adopt the given lower and upper bounds
  ! for the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd
  
  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  real(SP), dimension(:,:,:), pointer :: p_Sarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_SINGLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Fsingle3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_single3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Sarray => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_single3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double3DDef (ihandle, p_Darray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble3D

  end subroutine storage_getbase_double3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double3DUBnd (ihandle, p_Darray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array and adopt the given upper bound for the third dimension.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble3D,3))) then
       call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble3D(:,:,:ubnd)

  end subroutine storage_getbase_double3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_double3DLUBnd (ihandle, p_Darray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! double precision array and adopt the given lower and upper bounds
  ! for the third dimension.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(DP), dimension(:,:,:), pointer :: p_Darray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_DOUBLE) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Ddouble3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_double3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Darray => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_double3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quad3DDef (ihandle, p_Qarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad precision array.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(QP), dimension(:,:,:), pointer :: p_Qarray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad3D

  end subroutine storage_getbase_quad3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quad3DUBnd (ihandle, p_Qarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad precision array and adopt the given upper bound for the third dimension.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(QP), dimension(:,:,:), pointer :: p_Qarray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Qquad3D,3))) then
       call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad3D(:,:,:ubnd)

  end subroutine storage_getbase_quad3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_quad3DLUBnd (ihandle, p_Qarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! quad precision array and adopt the given lower and upper bounds
  ! for the third dimension.

!</description>

!<input>
  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!<output>
  ! The pointer associated to the handle.
  real(QP), dimension(:,:,:), pointer :: p_Qarray
!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_QUAD) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Qquad3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Qquad3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_quad3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Qarray => p_rheap%p_Rdescriptors(ihandle)%p_Qquad3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_quad3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical3DDef (ihandle, p_Larray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical3D

  end subroutine storage_getbase_logical3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical3DUBnd (ihandle, p_Larray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopt the given upper bound for the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Blogical3D,3))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical3D(:,:,:ubnd)

  end subroutine storage_getbase_logical3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_logical3DLUBnd (ihandle, p_Larray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! logical array and adopt the given lower and upper bounds for
  ! the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  logical, dimension(:,:,:), pointer :: p_Larray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_LOGICAL) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Blogical3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Blogical3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_logical3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Larray => p_rheap%p_Rdescriptors(ihandle)%p_Blogical3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_logical3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char3DDef (ihandle, p_Carray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar3D

  end subroutine storage_getbase_char3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char3DUBnd (ihandle, p_Carray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopt the given upper bound for the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Schar3D,3))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar3D(:,:,:ubnd)

  end subroutine storage_getbase_char3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_char3DLUBnd (ihandle, p_Carray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to a
  ! character array and adopt the given lower and upper bounds for
  ! the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  character, dimension(:,:,:), pointer :: p_Carray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_CHAR) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Schar3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Schar3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_char3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Carray => p_rheap%p_Rdescriptors(ihandle)%p_Schar3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_char3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8_3DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_3D

  end subroutine storage_getbase_int8_3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8_3DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_3D,3))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_3D(:,:,:ubnd)

  end subroutine storage_getbase_int8_3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int8_3DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I8), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT8) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint8_3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int8_3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_int8_3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16_3DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_3D

  end subroutine storage_getbase_int16_3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16_3DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the secondthird dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_3D,3))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_3D(:,:,:ubnd)

  end subroutine storage_getbase_int16_3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int16_3DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I16), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT16) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint16_3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int16_3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_int16_3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32_3DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_3D

  end subroutine storage_getbase_int32_3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32_3DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_3D,3))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_3D(:,:,:ubnd)

  end subroutine storage_getbase_int32_3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int32_3DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I32), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT32) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint32_3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_int32_3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64_3DDef (ihandle, p_Iarray, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_3DDef')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_3DDef')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_3D

  end subroutine storage_getbase_int64_3DDef

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64_3DUBnd (ihandle, p_Iarray, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given upper bound for the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_3DUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_3DUBnd')
      call sys_halt()
    end if

    if ((ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_3D,3))) then
      call output_line ('Upper bound invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_3D(:,:,:ubnd)

  end subroutine storage_getbase_int64_3DUBnd

!************************************************************************

!<subroutine>

  subroutine storage_getbase_int64_3DLUBnd (ihandle, p_Iarray, lbnd, ubnd, rheap)

!<description>

  ! This routine returns the pointer to a handle associated to an
  ! integer array and adopts the given lower and upper bounds for
  ! the third dimension.

!</description>

!<input>

  ! The handle
  integer, intent(in) :: ihandle

  ! The lower bound
  integer, intent(in) :: lbnd

  ! The upper bound
  integer, intent(in) :: ubnd

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

!</input>

!<output>

  ! The pointer associated to the handle.
  integer(I64), dimension(:,:,:), pointer :: p_Iarray

!</output>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_3DLUBnd')
      call sys_halt()
    end if

    if (p_rheap%p_Rdescriptors(ihandle)%idataType .ne. ST_INT64) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int32_3DLUBnd')
      call sys_halt()
    end if

    if ((lbnd .lt. 1) .or. &
        (lbnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_3D,3)) .or. &
        (ubnd .lt. 1) .or. &
        (ubnd .gt. ubound(p_rheap%p_Rdescriptors(ihandle)%p_Iint64_3D,3)) .or. &
        (lbnd .gt. ubnd)) then
      call output_line ('Bounds invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getbase_int64_3DUBnd')
      call sys_halt()
    end if

    ! Get the pointer
    p_Iarray => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_3D(:,:,lbnd:ubnd)

  end subroutine storage_getbase_int64_3DLUBnd

!************************************************************************

!<subroutine>

  subroutine storage_copyDefault (h_source, h_dest, rheap)

!<description>
  ! This routine copies the information of one array to another.
  ! The structure of the arrays behind h_source and h_dest must be the
  ! same!
!</description>

!<input>
  ! Handle of the source array to copy
  integer, intent(in) :: h_source

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  integer, intent(inout) :: h_dest
!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rsource, p_rdest
  integer :: isizeSource,isizeDest
  integer, dimension(2) :: Ilbound2D, Iubound2D
  integer, dimension(3) :: Ilbound3D, Iubound3D
  integer, dimension(2) :: Isize2DSource,Isize2DDest
  integer, dimension(3) :: Isize3DSource,Isize3DDest

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (h_source .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
      call sys_halt()
    end if
    if (.not. associated(p_rheap%p_Rdescriptors)) then
      call output_line ('Heap not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
      call sys_halt()
    end if

    p_rsource => p_rheap%p_Rdescriptors(h_source)
    
    ! Create a new array?
    if (h_dest .le. ST_NOHANDLE) then
      ! Create a new array in the same size and structure
      ! as h_source.
      select case (p_rsource%idimension)
      case (1)
        select case (p_rsource%idataType)
        case (ST_SINGLE)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Fsingle1D,1),&
                            ubound(p_rsource%p_Fsingle1D,1),&
                            ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_DOUBLE)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Ddouble1D,1),&
                            ubound(p_rsource%p_Ddouble1D,1),&
                            ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_QUAD)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Qquad1D,1),&
                            ubound(p_rsource%p_Qquad1D,1),&
                            ST_QUAD, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Iinteger1D,1),&
                            ubound(p_rsource%p_Iinteger1D,1),&
                            ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT8)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Iint8_1D,1),&
                            ubound(p_rsource%p_Iint8_1D,1),&
                            ST_INT8, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT16)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Iint16_1D,1),&
                            ubound(p_rsource%p_Iint16_1D,1),&
                            ST_INT16, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT32)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Iint32_1D,1),&
                            ubound(p_rsource%p_Iint32_1D,1),&
                            ST_INT32, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT64)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Iint64_1D,1),&
                            ubound(p_rsource%p_Iint64_1D,1),&
                            ST_INT64, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_LOGICAL)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Blogical1D,1),&
                            ubound(p_rsource%p_Blogical1D,1),&
                            ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_CHAR)
          call storage_new ('storage_copyDefault',p_rsource%sname,&
                            lbound(p_rsource%p_Schar1D,1),&
                            ubound(p_rsource%p_Schar1D,1),&
                            ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        end select
        
      case (2)
        select case (p_rsource%IdataType)
        case (ST_SINGLE)
          Ilbound2D = lbound(p_rsource%p_Fsingle2D)
          Iubound2D = ubound(p_rsource%p_Fsingle2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_DOUBLE)
          Ilbound2D = lbound(p_rsource%p_Ddouble2D)
          Iubound2D = ubound(p_rsource%p_Ddouble2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_QUAD)
          Ilbound2D = lbound(p_rsource%p_Qquad2D)
          Iubound2D = ubound(p_rsource%p_Qquad2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_QUAD, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT)
          Ilbound2D = lbound(p_rsource%p_Iinteger2D)
          Iubound2D = ubound(p_rsource%p_Iinteger2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT8)
          Ilbound2D = lbound(p_rsource%p_Iint8_2D)
          Iubound2D = ubound(p_rsource%p_Iint8_2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_INT8, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT16)
          Ilbound2D = lbound(p_rsource%p_Iint16_2D)
          Iubound2D = ubound(p_rsource%p_Iint16_2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_INT16, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT32)
          Ilbound2D = lbound(p_rsource%p_Iint32_2D)
          Iubound2D = ubound(p_rsource%p_Iint32_2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_INT32, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT64)
          Ilbound2D = lbound(p_rsource%p_Iint64_2D)
          Iubound2D = ubound(p_rsource%p_Iint64_2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_INT64, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_LOGICAL)
          Ilbound2D = lbound(p_rsource%p_Blogical2D)
          Iubound2D = ubound(p_rsource%p_Blogical2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_CHAR)
          Ilbound2D = lbound(p_rsource%p_Schar2D)
          Iubound2D = ubound(p_rsource%p_Schar2D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound2D,Iubound2D,&
                            ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        end select

      case (3)
        select case (p_rsource%IdataType)
        case (ST_SINGLE)
          Ilbound3D = lbound(p_rsource%p_Fsingle3D)
          Iubound3D = ubound(p_rsource%p_Fsingle3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_DOUBLE)
          Ilbound3D = lbound(p_rsource%p_Ddouble3D)
          Iubound3D = ubound(p_rsource%p_Ddouble3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_QUAD)
          Ilbound3D = lbound(p_rsource%p_Qquad3D)
          Iubound3D = ubound(p_rsource%p_Qquad3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_QUAD, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT)
          Ilbound3D = lbound(p_rsource%p_Iinteger3D)
          Iubound3D = ubound(p_rsource%p_Iinteger3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT8)
          Ilbound3D = lbound(p_rsource%p_Iint8_3D)
          Iubound3D = ubound(p_rsource%p_Iint8_3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_INT8, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT16)
          Ilbound3D = lbound(p_rsource%p_Iint16_3D)
          Iubound3D = ubound(p_rsource%p_Iint16_3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_INT16, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT32)
          Ilbound3D = lbound(p_rsource%p_Iint32_3D)
          Iubound3D = ubound(p_rsource%p_Iint32_3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_INT32, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_INT64)
          Ilbound3D = lbound(p_rsource%p_Iint64_3D)
          Iubound3D = ubound(p_rsource%p_Iint64_3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_INT64, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_LOGICAL)
          Ilbound3D = lbound(p_rsource%p_Blogical3D)
          Iubound3D = ubound(p_rsource%p_Blogical3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        case (ST_CHAR)
          Ilbound3D = lbound(p_rsource%p_Schar3D)
          Iubound3D = ubound(p_rsource%p_Schar3D)
          call storage_new ('storage_copyDefault', p_rsource%sname,&
                            Ilbound3D,Iubound3D,&
                            ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
        end select
      end select

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it is correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      p_rdest   => p_rheap%p_Rdescriptors(h_dest)

    else

      ! Check if the given destination handle is compatible with the source handle
      p_rdest => p_rheap%p_Rdescriptors(h_dest)

      ! Check if source and destination handle have the same dimension
      if (p_rsource%idimension .ne. p_rdest%idimension) then
        call output_line ('Dimension of source and destination handles are different!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
        call sys_halt()
      end if

      ! Sizes the same?
      select case(p_rsource%idimension)
      case (1)
        call storage_getsize(h_source, isizeSource, rheap)
        call storage_getsize(h_dest,   isizeDest,   rheap)
        if (isizeSource .ne. isizeDest) then
          call output_line ('Size of source and destination handles are different!', &
                             OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
        end if

      case (2)
        call storage_getsize(h_source, Isize2DSource, rheap)
        call storage_getsize(h_dest,   Isize2DDest,   rheap)
        if ((Isize2DSource(1) .ne. Isize2DDest(1)) .or.&
            (Isize2DSource(2) .ne. Isize2DDest(2))) then
          call output_line ('Size of source and destination handles are different!', &
                             OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
        end if

      case (3)
        call storage_getsize(h_source, Isize3DSource, rheap)
        call storage_getsize(h_dest,   Isize3DDest,   rheap)
        if ((Isize3DSource(1) .ne. Isize3DDest(1)) .or.&
            (Isize3DSource(2) .ne. Isize3DDest(2)) .or.&
            (Isize3DSource(3) .ne. Isize3DDest(3))) then
          call output_line ('Size of source and destination handles are different!', &
                             OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
        end if

      end select

    end if

    ! What is to copy
    select case (p_rsource%idimension)
    case (1)
      select case (p_rsource%idataType)
      case (ST_SINGLE)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorSngl (p_rsource%p_Fsingle1D,p_rdest%p_Fsingle1D)
        case (ST_DOUBLE)
          call lalg_copyVectorSnglDbl (p_rsource%p_Fsingle1D,p_rdest%p_Ddouble1D)
        case (ST_QUAD)
          call lalg_copyVectorSnglQuad (p_rsource%p_Fsingle1D,p_rdest%p_Qquad1D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_DOUBLE)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorDblSngl (p_rsource%p_Ddouble1D,p_rdest%p_Fsingle1D)
        case (ST_DOUBLE)
          call lalg_copyVectorDble (p_rsource%p_Ddouble1D,p_rdest%p_Ddouble1D)
        case (ST_QUAD)
          call lalg_copyVectorDblQuad (p_rsource%p_Ddouble1D,p_rdest%p_Qquad1D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_QUAD)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorQuadSngl (p_rsource%p_Qquad1D,p_rdest%p_Fsingle1D)
        case (ST_DOUBLE)
          call lalg_copyVectorQuadDbl (p_rsource%p_Qquad1D,p_rdest%p_Ddouble1D)
        case (ST_QUAD)
          call lalg_copyVectorQuad (p_rsource%p_Qquad1D,p_rdest%p_Qquad1D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt (p_rsource%p_Iinteger1D,p_rdest%p_Iinteger1D)
        case (ST_INT8)
          call lalg_copyVectorInt (p_rsource%p_Iinteger1D,p_rdest%p_Iint8_1D)
        case (ST_INT16)
          call lalg_copyVectorInt (p_rsource%p_Iinteger1D,p_rdest%p_Iint16_1D)
        case (ST_INT32)
          call lalg_copyVectorInt (p_rsource%p_Iinteger1D,p_rdest%p_Iint32_1D)
        case (ST_INT64)
          call lalg_copyVectorInt (p_rsource%p_Iinteger1D,p_rdest%p_Iint64_1D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT8)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt (p_rsource%p_Iint8_1D,p_rdest%p_Iinteger1D)
        case (ST_INT8)
          call lalg_copyVectorInt (p_rsource%p_Iint8_1D,p_rdest%p_Iint8_1D)
        case (ST_INT16)
          call lalg_copyVectorInt (p_rsource%p_Iint8_1D,p_rdest%p_Iint16_1D)
        case (ST_INT32)
          call lalg_copyVectorInt (p_rsource%p_Iint8_1D,p_rdest%p_Iint32_1D)
        case (ST_INT64)
          call lalg_copyVectorInt (p_rsource%p_Iint8_1D,p_rdest%p_Iint64_1D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT16)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt (p_rsource%p_Iint16_1D,p_rdest%p_Iinteger1D)
        case (ST_INT8)
          call lalg_copyVectorInt (p_rsource%p_Iint16_1D,p_rdest%p_Iint8_1D)
        case (ST_INT16)
          call lalg_copyVectorInt (p_rsource%p_Iint16_1D,p_rdest%p_Iint16_1D)
        case (ST_INT32)
          call lalg_copyVectorInt (p_rsource%p_Iint16_1D,p_rdest%p_Iint32_1D)
        case (ST_INT64)
          call lalg_copyVectorInt (p_rsource%p_Iint16_1D,p_rdest%p_Iint64_1D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT32)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt (p_rsource%p_Iint32_1D,p_rdest%p_Iinteger1D)
        case (ST_INT8)
          call lalg_copyVectorInt (p_rsource%p_Iint32_1D,p_rdest%p_Iint8_1D)
        case (ST_INT16)
          call lalg_copyVectorInt (p_rsource%p_Iint32_1D,p_rdest%p_Iint16_1D)
        case (ST_INT32)
          call lalg_copyVectorInt (p_rsource%p_Iint32_1D,p_rdest%p_Iint32_1D)
        case (ST_INT64)
          call lalg_copyVectorInt (p_rsource%p_Iint32_1D,p_rdest%p_Iint64_1D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT64)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt (p_rsource%p_Iint64_1D,p_rdest%p_Iinteger1D)
        case (ST_INT8)
          call lalg_copyVectorInt (p_rsource%p_Iint64_1D,p_rdest%p_Iint8_1D)
        case (ST_INT16)
          call lalg_copyVectorInt (p_rsource%p_Iint64_1D,p_rdest%p_Iint16_1D)
        case (ST_INT32)
          call lalg_copyVectorInt (p_rsource%p_Iint64_1D,p_rdest%p_Iint32_1D)
        case (ST_INT64)
          call lalg_copyVectorInt (p_rsource%p_Iint64_1D,p_rdest%p_Iint64_1D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_LOGICAL)
        if (p_rdest%idataType .eq. ST_LOGICAL) then
          call lalg_copyVectorLogical(p_rsource%p_Blogical1D,p_rdest%p_Blogical1D)
        else
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end if

      case (ST_CHAR)
        if (p_rdest%idataType .eq. ST_CHAR) then
          call lalg_copyVectorChar(p_rsource%p_Schar1D,p_rdest%p_Schar1D)
        else
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end if

      case default
        call output_line ('Unknown data type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
        call sys_halt()
      end select

    case (2)
      select case (p_rsource%idataType)
      case (ST_SINGLE)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorSngl2D (p_rsource%p_Fsingle2D,p_rdest%p_Fsingle2D)
        case (ST_DOUBLE)
          call lalg_copyVectorSnglDbl2D (p_rsource%p_Fsingle2D,p_rdest%p_Ddouble2D)
        case (ST_QUAD)
          call lalg_copyVectorSnglQuad2D (p_rsource%p_Fsingle2D,p_rdest%p_Qquad2D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select
        
      case (ST_DOUBLE)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorDblSngl2D (p_rsource%p_Ddouble2D,p_rdest%p_Fsingle2D)
        case (ST_DOUBLE)
          call lalg_copyVectorDble2D (p_rsource%p_Ddouble2D,p_rdest%p_Ddouble2D)
        case (ST_QUAD)
          call lalg_copyVectorDblQuad2D (p_rsource%p_Ddouble2D,p_rdest%p_Qquad2D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select
        
      case (ST_QUAD)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorQuadSngl2D (p_rsource%p_Qquad2D,p_rdest%p_Fsingle2D)
        case (ST_DOUBLE)
          call lalg_copyVectorQuadDbl2D (p_rsource%p_Qquad2D,p_rdest%p_Ddouble2D)
        case (ST_QUAD)
          call lalg_copyVectorQuad2D (p_rsource%p_Qquad2D,p_rdest%p_Qquad2D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt2D (p_rsource%p_Iinteger2D, p_rdest%p_Iinteger2D)
        case (ST_INT8)
          call lalg_copyVectorInt2D (p_rsource%p_Iinteger2D, p_rdest%p_Iint8_2D)
        case (ST_INT16)
          call lalg_copyVectorInt2D (p_rsource%p_Iinteger2D, p_rdest%p_Iint16_2D)
        case (ST_INT32)
          call lalg_copyVectorInt2D (p_rsource%p_Iinteger2D, p_rdest%p_Iint32_2D)
        case (ST_INT64)
          call lalg_copyVectorInt2D (p_rsource%p_Iinteger2D, p_rdest%p_Iint64_2D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT8)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt2D (p_rsource%p_Iint8_2D, p_rdest%p_Iinteger2D)
        case (ST_INT8)
          call lalg_copyVectorInt2D (p_rsource%p_Iint8_2D, p_rdest%p_Iint8_2D)
        case (ST_INT16)
          call lalg_copyVectorInt2D (p_rsource%p_Iint8_2D, p_rdest%p_Iint16_2D)
        case (ST_INT32)
          call lalg_copyVectorInt2D (p_rsource%p_Iint8_2D, p_rdest%p_Iint32_2D)
        case (ST_INT64)
          call lalg_copyVectorInt2D (p_rsource%p_Iint8_2D, p_rdest%p_Iint64_2D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT16)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt2D (p_rsource%p_Iint16_2D, p_rdest%p_Iinteger2D)
        case (ST_INT8)
          call lalg_copyVectorInt2D (p_rsource%p_Iint16_2D, p_rdest%p_Iint8_2D)
        case (ST_INT16)
          call lalg_copyVectorInt2D (p_rsource%p_Iint16_2D, p_rdest%p_Iint16_2D)
        case (ST_INT32)
          call lalg_copyVectorInt2D (p_rsource%p_Iint16_2D, p_rdest%p_Iint32_2D)
        case (ST_INT64)
          call lalg_copyVectorInt2D (p_rsource%p_Iint16_2D, p_rdest%p_Iint64_2D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT32)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt2D (p_rsource%p_Iint32_2D, p_rdest%p_Iinteger2D)
        case (ST_INT8)
          call lalg_copyVectorInt2D (p_rsource%p_Iint32_2D, p_rdest%p_Iint8_2D)
        case (ST_INT16)
          call lalg_copyVectorInt2D (p_rsource%p_Iint32_2D, p_rdest%p_Iint16_2D)
        case (ST_INT32)
          call lalg_copyVectorInt2D (p_rsource%p_Iint32_2D, p_rdest%p_Iint32_2D)
        case (ST_INT64)
          call lalg_copyVectorInt2D (p_rsource%p_Iint32_2D, p_rdest%p_Iint64_2D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT64)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt2D (p_rsource%p_Iint64_2D, p_rdest%p_Iinteger2D)
        case (ST_INT8)
          call lalg_copyVectorInt2D (p_rsource%p_Iint64_2D, p_rdest%p_Iint8_2D)
        case (ST_INT16)
          call lalg_copyVectorInt2D (p_rsource%p_Iint64_2D, p_rdest%p_Iint16_2D)
        case (ST_INT32)
          call lalg_copyVectorInt2D (p_rsource%p_Iint64_2D, p_rdest%p_Iint32_2D)
        case (ST_INT64)
          call lalg_copyVectorInt2D (p_rsource%p_Iint64_2D, p_rdest%p_Iint64_2D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_LOGICAL)
        if (p_rdest%idataType .eq. ST_LOGICAL) then
          call lalg_copyVectorLogical2D (p_rsource%p_Blogical2D, p_rdest%p_Blogical2D)
        else
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end if

      case (ST_CHAR)
        if (p_rdest%idataType .eq. ST_CHAR) then
          call lalg_copyVectorChar2D (p_rsource%p_Schar2D, p_rdest%p_Schar2D)
        else
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end if

      case default
        call output_line ('Unknown data type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
        call sys_halt()
      end select

    case (3)
      select case (p_rsource%idataType)
      case (ST_SINGLE)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorSngl3D (p_rsource%p_Fsingle3D,p_rdest%p_Fsingle3D)
        case (ST_DOUBLE)
          call lalg_copyVectorSnglDbl3D (p_rsource%p_Fsingle3D,p_rdest%p_Ddouble3D)
        case (ST_QUAD)
          call lalg_copyVectorSnglQuad3D (p_rsource%p_Fsingle3D,p_rdest%p_Qquad3D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select
        
      case (ST_DOUBLE)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorDblSngl3D (p_rsource%p_Ddouble3D,p_rdest%p_Fsingle3D)
        case (ST_DOUBLE)
          call lalg_copyVectorDble3D (p_rsource%p_Ddouble3D,p_rdest%p_Ddouble3D)
        case (ST_QUAD)
          call lalg_copyVectorDblQuad3D (p_rsource%p_Ddouble3D,p_rdest%p_Qquad3D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select
        
      case (ST_QUAD)
        select case (p_rdest%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorQuadSngl3D (p_rsource%p_Qquad3D,p_rdest%p_Fsingle3D)
        case (ST_DOUBLE)
          call lalg_copyVectorQuadDbl3D (p_rsource%p_Qquad3D,p_rdest%p_Ddouble3D)
        case (ST_QUAD)
          call lalg_copyVectorQuad3D (p_rsource%p_Qquad3D,p_rdest%p_Qquad3D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt3D (p_rsource%p_Iinteger3D, p_rdest%p_Iinteger3D)
        case (ST_INT8)
          call lalg_copyVectorInt3D (p_rsource%p_Iinteger3D, p_rdest%p_Iint8_3D)
        case (ST_INT16)
          call lalg_copyVectorInt3D (p_rsource%p_Iinteger3D, p_rdest%p_Iint16_3D)
        case (ST_INT32)
          call lalg_copyVectorInt3D (p_rsource%p_Iinteger3D, p_rdest%p_Iint32_3D)
        case (ST_INT64)
          call lalg_copyVectorInt3D (p_rsource%p_Iinteger3D, p_rdest%p_Iint64_3D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT8)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt3D (p_rsource%p_Iint8_3D, p_rdest%p_Iinteger3D)
        case (ST_INT8)
          call lalg_copyVectorInt3D (p_rsource%p_Iint8_3D, p_rdest%p_Iint8_3D)
        case (ST_INT16)
          call lalg_copyVectorInt3D (p_rsource%p_Iint8_3D, p_rdest%p_Iint16_3D)
        case (ST_INT32)
          call lalg_copyVectorInt3D (p_rsource%p_Iint8_3D, p_rdest%p_Iint32_3D)
        case (ST_INT64)
          call lalg_copyVectorInt3D (p_rsource%p_Iint8_3D, p_rdest%p_Iint64_3D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT16)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt3D (p_rsource%p_Iint16_3D, p_rdest%p_Iinteger3D)
        case (ST_INT8)
          call lalg_copyVectorInt3D (p_rsource%p_Iint16_3D, p_rdest%p_Iint8_3D)
        case (ST_INT16)
          call lalg_copyVectorInt3D (p_rsource%p_Iint16_3D, p_rdest%p_Iint16_3D)
        case (ST_INT32)
          call lalg_copyVectorInt3D (p_rsource%p_Iint16_3D, p_rdest%p_Iint32_3D)
        case (ST_INT64)
          call lalg_copyVectorInt3D (p_rsource%p_Iint16_3D, p_rdest%p_Iint64_3D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT32)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt3D (p_rsource%p_Iint32_3D, p_rdest%p_Iinteger3D)
        case (ST_INT8)
          call lalg_copyVectorInt3D (p_rsource%p_Iint32_3D, p_rdest%p_Iint8_3D)
        case (ST_INT16)
          call lalg_copyVectorInt3D (p_rsource%p_Iint32_3D, p_rdest%p_Iint16_3D)
        case (ST_INT32)
          call lalg_copyVectorInt3D (p_rsource%p_Iint32_3D, p_rdest%p_Iint32_3D)
        case (ST_INT64)
          call lalg_copyVectorInt3D (p_rsource%p_Iint32_3D, p_rdest%p_Iint64_3D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_INT64)
        select case (p_rdest%idataType)
        case (ST_INT)
          call lalg_copyVectorInt3D (p_rsource%p_Iint64_3D, p_rdest%p_Iinteger3D)
        case (ST_INT8)
          call lalg_copyVectorInt3D (p_rsource%p_Iint64_3D, p_rdest%p_Iint8_3D)
        case (ST_INT16)
          call lalg_copyVectorInt3D (p_rsource%p_Iint64_3D, p_rdest%p_Iint16_3D)
        case (ST_INT32)
          call lalg_copyVectorInt3D (p_rsource%p_Iint64_3D, p_rdest%p_Iint32_3D)
        case (ST_INT64)
          call lalg_copyVectorInt3D (p_rsource%p_Iint64_3D, p_rdest%p_Iint64_3D)
        case default
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end select

      case (ST_LOGICAL)
        if (p_rdest%idataType .eq. ST_LOGICAL) then
          call lalg_copyVectorLogical3D (p_rsource%p_Blogical3D, p_rdest%p_Blogical3D)
        else
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end if

      case (ST_CHAR)
        if (p_rdest%idataType .eq. ST_CHAR) then
          call lalg_copyVectorChar3D (p_rsource%p_Schar3D, p_rdest%p_Schar3D)
        else
          call output_line ('Unsupported data type combination!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
          call sys_halt()
        end if

      case default
        call output_line ('Unknown data type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copyDefault')
        call sys_halt()
      end select
      
    end select

  end subroutine storage_copyDefault

!************************************************************************

!<subroutine>

  subroutine storage_copy_explicit (h_source, h_dest, Istart_source, &
                                    Istart_dest, Ilength, rheap)

!<description>
  ! This routine copies the information of one array to another.
  ! The structure of the arrays behind h_source and h_dest must be
  ! "similar" in the following sense: Both datatypes and dimensions
  ! must be the same. The routine allows for copying only parts of
  ! of the arrays. Therefor the relevant parts of the arrays must
  ! be the same!
  ! It is actually a wrapper for subroutines 'storage_copy_explicitx"'.
!</description>

!<input>
  ! Handle of the source array to copy
  integer, intent(in) :: h_source
  
  ! First entry of the source array to copy
  integer, dimension(:), intent(in) :: Istart_source

  ! First entry of the destination array where to copy
  integer, dimension(:), intent(in) :: Istart_dest

  ! Length of the array to copy
  integer, dimension(:), intent(in) :: Ilength

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  integer, intent(inout) :: h_dest
!</inputoutput>

!</subroutine>

    if (size(Istart_source) .ne. size(Istart_dest) .or.&
        size(Istart_source) .ne. size(Ilength)) then
      call output_line ('dim(Istart_source) /= dim(Istart_dest) /= dim(Ilength)!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit')
      call sys_halt()
    end if
    
    select case (size(Ilength))
    case (1)
      call storage_copy_explicit1d (h_source, h_dest, Istart_source(1),&
                                    Istart_dest(1), Ilength(1), rheap)
    case (2)
      call storage_copy_explicit2d (h_source, h_dest, Istart_source,&
                                    Istart_dest, Ilength, rheap)
    case (3)
      call storage_copy_explicit3d (h_source, h_dest, Istart_source,&
                                    Istart_dest, Ilength, rheap)
    case default
      call output_line ('Memory blocks of dimensions larger than 3 are not supported', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit')
    end select

  end subroutine storage_copy_explicit

!************************************************************************

!<subroutine>

  subroutine storage_copy_explicit1d (h_source, h_dest, istart_source, &
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
  integer, intent(in) :: h_source

  ! First entry of the source array to copy
  integer, intent(in) :: istart_source

  ! First entry of the destination array where to copy
  integer, intent(in) :: istart_dest

  ! Length of the array to copy
  integer, intent(in) :: ilength

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  integer, intent(inout) :: h_dest
!</inputoutput>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rsource, p_rdest
  integer :: i

!!$    ! Check if the start address is positive
!!$    if (istart_source .le. 0 .or. istart_dest .le. 0) then
!!$      call output_line ('Start address must be positive!', &
!!$                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
!!$      call sys_halt()
!!$    end if

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (h_source .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
      call sys_halt()
    end if
    if (.not. associated(p_rheap%p_Rdescriptors)) then
      call output_line ('Heap not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
      call sys_halt()
    end if

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    if (h_dest .le. ST_NOHANDLE) then
      ! Create a new array in the same size and structure
      ! as h_source.
      if (p_rsource%idimension .ne. 1) then
        call output_line ('Only 1D arrays are allowed!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end if

      select case (p_rsource%idataType)
      case (ST_SINGLE)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Fsingle1D,1),&
                          ubound(p_rsource%p_Fsingle1D,1),&
                          ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_DOUBLE)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Ddouble1D,1),&
                          ubound(p_rsource%p_Ddouble1D,1),&
                          ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_QUAD)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Qquad1D,1),&
                          ubound(p_rsource%p_Qquad1D,1),&
                          ST_QUAD, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Iinteger1D,1),&
                          ubound(p_rsource%p_Iinteger1D,1),&
                          ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT8)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Iint8_1D,1),&
                          ubound(p_rsource%p_Iint8_1D,1),&
                          ST_INT8, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT16)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Iint16_1D,1),&
                          ubound(p_rsource%p_Iint16_1D,1),&
                          ST_INT16, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT32)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Iint32_1D,1),&
                          ubound(p_rsource%p_Iint32_1D,1),&
                          ST_INT32, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT64)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Iint64_1D,1),&
                          ubound(p_rsource%p_Iint64_1D,1),&
                          ST_INT64, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_LOGICAL)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Blogical1D,1),&
                          ubound(p_rsource%p_Blogical1D,1),&
                          ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_CHAR)
        call storage_new ('storage_copy_explicit1d',p_rsource%sname,&
                          lbound(p_rsource%p_Schar1D,1),&
                          ubound(p_rsource%p_Schar1D,1),&
                          ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      end select

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it is correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      p_rdest   => p_rheap%p_Rdescriptors(h_dest)

    else
      
      ! Check if the given destination handle is compatible with the source handle
      p_rdest => p_rheap%p_Rdescriptors(h_dest)

      ! 1D/2D the same?
      if (p_rsource%idimension .ne. p_rdest%idimension) then
        call output_line ('Dimension of source and destination handles are different!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end if

    end if
    
    ! What is to copy
    select case (p_rsource%idataType)
    case (ST_SINGLE)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (istart_source           < lbound(p_rsource%p_Fsingle1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Fsingle1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Fsingle1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Fsingle1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Fsingle1D(istart_dest+i-1) = &
              p_rsource%p_Fsingle1D(istart_source+i-1)
        end do

      case (ST_DOUBLE)
        if (istart_source           < lbound(p_rsource%p_Fsingle1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Ddouble1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Fsingle1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Ddouble1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Ddouble1D(istart_dest+i-1) = &
              p_rsource%p_Fsingle1D(istart_source+i-1)
        end do

      case (ST_QUAD)
        if (istart_source           < lbound(p_rsource%p_Fsingle1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Qquad1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Fsingle1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Qquad1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Qquad1D(istart_dest+i-1) = &
              p_rsource%p_Fsingle1D(istart_source+i-1)
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end select
      
    case (ST_DOUBLE)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (istart_source           < lbound(p_rsource%p_Ddouble1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Fsingle1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Ddouble1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Fsingle1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Fsingle1D(istart_dest+i-1) = &
              p_rsource%p_Ddouble1D(istart_source+i-1)
        end do

      case (ST_DOUBLE)
        if (istart_source           < lbound(p_rsource%p_Ddouble1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Ddouble1D,1) .or. &
            istart_source+ilength-1 > ubound(p_rsource%p_Ddouble1D,1) .or. &
            istart_dest+ilength-1   > ubound(p_rdest%p_Ddouble1D,1)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Ddouble1D(istart_dest+i-1) = &
              p_rsource%p_Ddouble1D(istart_source+i-1)
        end do

      case (ST_QUAD)
        if (istart_source           < lbound(p_rsource%p_Ddouble1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Qquad1D,1) .or. &
            istart_source+ilength-1 > ubound(p_rsource%p_Ddouble1D,1) .or. &
            istart_dest+ilength-1   > ubound(p_rdest%p_Qquad1D,1)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Qquad1D(istart_dest+i-1) = &
              p_rsource%p_Ddouble1D(istart_source+i-1)
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end select
      
    case (ST_QUAD)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (istart_source           < lbound(p_rsource%p_Qquad1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Fsingle1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Qquad1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Fsingle1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Fsingle1D(istart_dest+i-1) = &
              p_rsource%p_Qquad1D(istart_source+i-1)
        end do

      case (ST_DOUBLE)
        if (istart_source           < lbound(p_rsource%p_Qquad1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Ddouble1D,1) .or. &
            istart_source+ilength-1 > ubound(p_rsource%p_Qquad1D,1) .or. &
            istart_dest+ilength-1   > ubound(p_rdest%p_Ddouble1D,1)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Ddouble1D(istart_dest+i-1) = &
              p_rsource%p_Qquad1D(istart_source+i-1)
        end do

      case (ST_QUAD)
        if (istart_source           < lbound(p_rsource%p_Qquad1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Qquad1D,1) .or. &
            istart_source+ilength-1 > ubound(p_rsource%p_Qquad1D,1) .or. &
            istart_dest+ilength-1   > ubound(p_rdest%p_Qquad1D,1)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Qquad1D(istart_dest+i-1) = &
              p_rsource%p_Qquad1D(istart_source+i-1)
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end select
      
    case (ST_INT)
      select case (p_rdest%idataType)
      case (ST_INT)
        if (istart_source           < lbound(p_rsource%p_Iinteger1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iinteger1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iinteger1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iinteger1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iinteger1D(istart_dest+i-1) = &
              p_rsource%p_Iinteger1D(istart_source+i-1)
        end do

      case (ST_INT8)
        if (istart_source           < lbound(p_rsource%p_Iinteger1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint8_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iinteger1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint8_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint8_1D(istart_dest+i-1) = &
              p_rsource%p_Iinteger1D(istart_source+i-1)
        end do

      case (ST_INT16)
        if (istart_source           < lbound(p_rsource%p_Iinteger1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint16_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iinteger1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint16_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint16_1D(istart_dest+i-1) = &
              p_rsource%p_Iinteger1D(istart_source+i-1)
        end do

      case (ST_INT32)
        if (istart_source           < lbound(p_rsource%p_Iinteger1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint32_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iinteger1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint32_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint32_1D(istart_dest+i-1) = &
              p_rsource%p_Iinteger1D(istart_source+i-1)
        end do

      case (ST_INT64)
        if (istart_source           < lbound(p_rsource%p_Iinteger1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint64_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iinteger1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint64_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint64_1D(istart_dest+i-1) = &
              p_rsource%p_Iinteger1D(istart_source+i-1)
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end select
      
    case (ST_INT8)
      select case (p_rdest%idataType)
      case (ST_INT)
        if (istart_source           < lbound(p_rsource%p_Iint8_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iinteger1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint8_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iinteger1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iinteger1D(istart_dest+i-1) = &
              p_rsource%p_Iint8_1D(istart_source+i-1)
        end do

      case (ST_INT8)
        if (istart_source           < lbound(p_rsource%p_Iint8_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint8_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint8_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint8_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint8_1D(istart_dest+i-1) = &
              p_rsource%p_Iint8_1D(istart_source+i-1)
        end do

      case (ST_INT16)
        if (istart_source           < lbound(p_rsource%p_Iint8_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint16_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint8_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint16_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint16_1D(istart_dest+i-1) = &
              p_rsource%p_Iint8_1D(istart_source+i-1)
        end do

      case (ST_INT32)
        if (istart_source           < lbound(p_rsource%p_Iint8_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint32_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint8_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint32_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint32_1D(istart_dest+i-1) = &
              p_rsource%p_Iint8_1D(istart_source+i-1)
        end do

      case (ST_INT64)
        if (istart_source           < lbound(p_rsource%p_Iint8_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint64_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint8_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint64_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint64_1D(istart_dest+i-1) = &
              p_rsource%p_Iint8_1D(istart_source+i-1)
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end select
      
    case (ST_INT16)
      select case (p_rdest%idataType)
      case (ST_INT)
        if (istart_source           < lbound(p_rsource%p_Iint16_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iinteger1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint16_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iinteger1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iinteger1D(istart_dest+i-1) = &
              p_rsource%p_Iint16_1D(istart_source+i-1)
        end do

      case (ST_INT8)
        if (istart_source           < lbound(p_rsource%p_Iint16_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint8_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint16_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint8_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint8_1D(istart_dest+i-1) = &
              p_rsource%p_Iint16_1D(istart_source+i-1)
        end do

      case (ST_INT16)
        if (istart_source           < lbound(p_rsource%p_Iint16_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint16_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint16_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint16_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint16_1D(istart_dest+i-1) = &
              p_rsource%p_Iint16_1D(istart_source+i-1)
        end do

      case (ST_INT32)
        if (istart_source           < lbound(p_rsource%p_Iint16_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint32_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint16_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint32_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint32_1D(istart_dest+i-1) = &
              p_rsource%p_Iint16_1D(istart_source+i-1)
        end do

      case (ST_INT64)
        if (istart_source           < lbound(p_rsource%p_Iint16_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint64_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint16_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint64_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint64_1D(istart_dest+i-1) = &
              p_rsource%p_Iint16_1D(istart_source+i-1)
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end select
      
    case (ST_INT32)
      select case (p_rdest%idataType)
      case (ST_INT)
        if (istart_source           < lbound(p_rsource%p_Iint32_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iinteger1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint32_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iinteger1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iinteger1D(istart_dest+i-1) = &
              p_rsource%p_Iint32_1D(istart_source+i-1)
        end do

      case (ST_INT8)
        if (istart_source           < lbound(p_rsource%p_Iint32_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint8_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint32_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint8_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint8_1D(istart_dest+i-1) = &
              p_rsource%p_Iint32_1D(istart_source+i-1)
        end do

      case (ST_INT16)
        if (istart_source           < lbound(p_rsource%p_Iint32_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint16_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint32_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint16_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint16_1D(istart_dest+i-1) = &
              p_rsource%p_Iint32_1D(istart_source+i-1)
        end do

      case (ST_INT32)
        if (istart_source           < lbound(p_rsource%p_Iint32_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint32_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint32_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint32_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint32_1D(istart_dest+i-1) = &
              p_rsource%p_Iint32_1D(istart_source+i-1)
        end do

      case (ST_INT64)
        if (istart_source           < lbound(p_rsource%p_Iint32_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint64_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint32_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint64_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint64_1D(istart_dest+i-1) = &
              p_rsource%p_Iint32_1D(istart_source+i-1)
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end select
      
    case (ST_INT64)
      select case (p_rdest%idataType)
      case (ST_INT)
        if (istart_source           < lbound(p_rsource%p_Iint64_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iinteger1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint64_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iinteger1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iinteger1D(istart_dest+i-1) = &
              p_rsource%p_Iint64_1D(istart_source+i-1)
        end do

      case (ST_INT8)
        if (istart_source           < lbound(p_rsource%p_Iint64_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint8_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint64_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint8_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint8_1D(istart_dest+i-1) = &
              p_rsource%p_Iint64_1D(istart_source+i-1)
        end do

      case (ST_INT16)
        if (istart_source           < lbound(p_rsource%p_Iint64_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint16_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint64_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint16_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint16_1D(istart_dest+i-1) = &
              p_rsource%p_Iint64_1D(istart_source+i-1)
        end do

      case (ST_INT32)
        if (istart_source           < lbound(p_rsource%p_Iint64_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint32_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint64_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint32_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint32_1D(istart_dest+i-1) = &
              p_rsource%p_Iint64_1D(istart_source+i-1)
        end do

      case (ST_INT64)
        if (istart_source           < lbound(p_rsource%p_Iint64_1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Iint64_1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Iint64_1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Iint64_1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Iint64_1D(istart_dest+i-1) = &
              p_rsource%p_Iint64_1D(istart_source+i-1)
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end select

    case (ST_LOGICAL)
      if (p_rdest%idataType .eq. ST_LOGICAL) then
        if (istart_source           < lbound(p_rsource%p_Blogical1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Blogical1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Blogical1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Blogical1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Blogical1D(istart_dest+i-1) = &
              p_rsource%p_Blogical1D(istart_source+i-1)
        end do
      else
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end if
      
    case (ST_CHAR)
      if (p_rdest%idataType .eq. ST_CHAR) then
        if (istart_source           < lbound(p_rsource%p_Schar1D,1) .or. &
            istart_dest             < lbound(p_rdest%p_Schar1D,1) .or. &
            istart_source+ilength-1 > size(p_rsource%p_Schar1D) .or. &
            istart_dest+ilength-1   > size(p_rdest%p_Schar1D)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
          call sys_halt()
        end if
        ! Copy by hand
        do i=1,ilength
          p_rdest%p_Schar1D(istart_dest+i-1) = &
              p_rsource%p_Schar1D(istart_source+i-1)
        end do
      else
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
        call sys_halt()
      end if
      
    case default
      call output_line ('Unknown data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit1d')
      call sys_halt()
    end select
    
  end subroutine storage_copy_explicit1d

!************************************************************************

!<subroutine>

  subroutine storage_copy_explicit2D (h_source, h_dest, Istart_source, &
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
  integer, intent(in) :: h_source

  ! First entry of the source array to copy
  integer, dimension(2), intent(in) :: Istart_source

  ! First entry of the destination array where to copy
  integer, dimension(2), intent(in) :: Istart_dest

  ! Length of the array to copy
  integer, dimension(2), intent(in) :: Ilength

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  integer, intent(inout) :: h_dest
!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rsource, p_rdest
  integer :: i,j

!!$    ! Check if the start address is positive
!!$    if (any(istart_source .le. 0) .or. any(istart_dest .le. 0)) then
!!$      call output_line ('Start address must be positive!', &
!!$                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
!!$      call sys_halt()
!!$    end if

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (h_source .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
      call sys_halt()
    end if
    if (.not. associated(p_rheap%p_Rdescriptors)) then
      call output_line ('Heap not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
      call sys_halt()
    end if

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    if (h_dest .le. ST_NOHANDLE) then
      ! Create a new array in the same size and structure
      ! as h_source.
      if (p_rsource%idimension .ne. 2) then
        call output_line ('Only 2D arrays are allowed!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
      end if

      select case (p_rsource%IdataType)
      case (ST_SINGLE)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Fsingle2D),&
                          ubound(p_rsource%p_Fsingle2D),&
                          ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_DOUBLE)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Ddouble2D),&
                          ubound(p_rsource%p_Ddouble2D),&
                          ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_QUAD)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Qquad2D),&
                          ubound(p_rsource%p_Qquad2D),&
                          ST_QUAD, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Iinteger2D),&
                          ubound(p_rsource%p_Iinteger2D),&
                          ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT8)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Iint8_2D),&
                          ubound(p_rsource%p_Iint8_2D),&
                          ST_INT8, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT16)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Iint16_2D),&
                          ubound(p_rsource%p_Iint16_2D),&
                          ST_INT16, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT32)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Iint32_2D),&
                          ubound(p_rsource%p_Iint32_2D),&
                          ST_INT32, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT64)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Iint64_2D),&
                          ubound(p_rsource%p_Iint64_2D),&
                          ST_INT64, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_LOGICAL)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Blogical2D),&
                          ubound(p_rsource%p_Blogical2D),&
                          ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_CHAR)
        call storage_new ('storage_copy_explicit2D', p_rsource%sname,&
                          lbound(p_rsource%p_Schar2D),&
                          ubound(p_rsource%p_Schar2D),&
                          ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      end select

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it is correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      p_rdest => p_rheap%p_Rdescriptors(h_dest)

    else
      
      ! Check if the given destination handle is compatible with the source handle
      p_rdest => p_rheap%p_Rdescriptors(h_dest)
      
      ! 1D/2D the same?
      if (p_rsource%idimension .ne. p_rdest%idimension) then
        call output_line ('Dimension of source and destination handles are different!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
      end if

    end if

    ! What is to copy
    select case (p_rsource%idataType)
    case (ST_SINGLE)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Fsingle2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Fsingle2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Fsingle2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Fsingle2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Fsingle2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_DOUBLE)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Ddouble2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Ddouble2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Ddouble2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Ddouble2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Ddouble2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_QUAD)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Qquad2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Qquad2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Qquad2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Qquad2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Qquad2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint8_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint16_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint32_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint64_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Fsingle2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
        
      end select
      
    case (ST_DOUBLE)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Fsingle2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Fsingle2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Fsingle2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Fsingle2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Fsingle2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_DOUBLE)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Ddouble2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Ddouble2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Ddouble2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Ddouble2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Ddouble2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_QUAD)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Qquad2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Qquad2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Qquad2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Qquad2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Qquad2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint8_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint16_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint32_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do

      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint64_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Ddouble2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
        
      end select
      
    case (ST_QUAD)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Fsingle2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Fsingle2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Fsingle2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Fsingle2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Fsingle2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Qquad2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_DOUBLE)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Ddouble2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Ddouble2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Ddouble2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Ddouble2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Ddouble2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Qquad2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_QUAD)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Qquad2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Qquad2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Qquad2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Qquad2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Qquad2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Qquad2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Qquad2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint8_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Qquad2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint16_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Qquad2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint32_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Qquad2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do

      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint64_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Qquad2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
        
      end select
      
    case (ST_INT)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iinteger2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint8_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iinteger2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint16_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iinteger2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint32_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iinteger2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint64_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iinteger2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
        
      end select
      
    case (ST_INT8)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint8_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint8_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint8_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint16_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint8_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint32_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint8_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint64_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint8_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
        
      end select
      
    case (ST_INT16)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint16_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint8_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint16_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint16_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint16_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint32_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint16_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint64_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint16_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
        
      end select
      
    case (ST_INT32)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint32_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint8_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint32_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint16_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint32_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint32_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint32_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint64_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint32_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
        
      end select
      
    case (ST_INT64)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iinteger2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint64_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint8_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint64_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint16_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint64_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint32_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint64_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_2D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_2D,2) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_2D,2) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_2D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_2D,2) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_2D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_2D,2)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
          call sys_halt()
        end if
        ! Copy by hand
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Iint64_2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
                p_rsource%p_Iint64_2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
        
      end select


    case (ST_LOGICAL)
      if (p_rdest%idataType .ne. ST_LOGICAL) then
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
      end if
      if (Istart_source(1) < lbound(p_rsource%p_Blogical2D,1) .or.&
          Istart_source(2) < lbound(p_rsource%p_Blogical2D,2) .or.&
          Istart_dest(1)   < lbound(p_rdest%p_Blogical2D,1) .or. &
          Istart_dest(2)   < lbound(p_rdest%p_Blogical2D,2) .or. &
          Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Blogical2D,1) .or. &
          Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Blogical2D,2) .or. &
          Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Blogical2D,1) .or. &
          Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Blogical2D,2)) then
        call output_line ('Subarrays incompatible!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
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
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
      end if
      if (Istart_source(1) < lbound(p_rsource%p_Schar2D,1) .or.&
          Istart_source(2) < lbound(p_rsource%p_Schar2D,2) .or.&
          Istart_dest(1)   < lbound(p_rdest%p_Schar2D,1) .or. &
          Istart_dest(2)   < lbound(p_rdest%p_Schar2D,2) .or. &
          Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Schar2D,1) .or. &
          Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Schar2D,2) .or. &
          Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Schar2D,1) .or. &
          Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Schar2D,2)) then
        call output_line ('Subarrays incompatible!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
        call sys_halt()
      end if
      ! Copy by hand
      do j=1,Ilength(2)
        do i=1,Ilength(1)
          p_rdest%p_Schar2D(Istart_dest(1)+i-1,Istart_dest(2)+j-1) = &
              p_rsource%p_Schar2D(Istart_source(1)+i-1,Istart_source(2)+j-1)
        end do
      end do
      
    case default
      call output_line ('Unknown data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit2D')
      call sys_halt()
    end select

  end subroutine storage_copy_explicit2D

!************************************************************************

!<subroutine>

  subroutine storage_copy_explicit3D (h_source, h_dest, Istart_source, &
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
  integer, intent(in) :: h_source

  ! First entry of the source array to copy
  integer, dimension(3), intent(in) :: Istart_source

  ! First entry of the destination array where to copy
  integer, dimension(3), intent(in) :: Istart_dest

  ! Length of the array to copy
  integer, dimension(3), intent(in) :: Ilength

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<inputoutput>
  ! Handle of the destination array.
  ! If =ST_NOHANDLE, a new handle is allocated in exactly the same size
  ! and structure as h_source and data is copied to it.
  integer, intent(inout) :: h_dest
!</inputoutput>

!</subroutine>

  ! local variables

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rsource, p_rdest
  integer :: i,j,k

!!$    ! Check if the start address is positive
!!$    if (any(istart_source .le. 0) .or. any(istart_dest .le. 0)) then
!!$      call output_line ('Start address must be positive!', &
!!$                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
!!$      call sys_halt()
!!$    end if

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (h_source .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
      call sys_halt()
    end if
    if (.not. associated(p_rheap%p_Rdescriptors)) then
      call output_line ('Heap not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
      call sys_halt()
    end if

    p_rsource => p_rheap%p_Rdescriptors(h_source)

    ! Create a new array?
    if (h_dest .le. ST_NOHANDLE) then
      ! Create a new array in the same size and structure
      ! as h_source.
      if (p_rsource%idimension .ne. 3) then
        call output_line ('Only 3D arrays are allowed!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
      end if

      select case (p_rsource%IdataType)
      case (ST_SINGLE)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Fsingle3D),&
                          ubound(p_rsource%p_Fsingle3D),&
                          ST_SINGLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_DOUBLE)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Ddouble3D),&
                          ubound(p_rsource%p_Ddouble3D),&
                          ST_DOUBLE, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_QUAD)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Qquad3D),&
                          ubound(p_rsource%p_Qquad3D),&
                          ST_QUAD, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Iinteger3D),&
                          ubound(p_rsource%p_Iinteger3D),&
                          ST_INT, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT8)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Iint8_3D),&
                          ubound(p_rsource%p_Iint8_3D),&
                          ST_INT8, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT16)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Iint16_3D),&
                          ubound(p_rsource%p_Iint16_3D),&
                          ST_INT16, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT32)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Iint32_3D),&
                          ubound(p_rsource%p_Iint32_3D),&
                          ST_INT32, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_INT64)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Iint64_3D),&
                          ubound(p_rsource%p_Iint64_3D),&
                          ST_INT64, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_LOGICAL)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Blogical3D),&
                          ubound(p_rsource%p_Blogical3D),&
                          ST_LOGICAL, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      case (ST_CHAR)
        call storage_new ('storage_copy_explicit3D', p_rsource%sname,&
                          lbound(p_rsource%p_Schar3D),&
                          ubound(p_rsource%p_Schar3D),&
                          ST_CHAR, h_dest, ST_NEWBLOCK_NOINIT, p_rheap)
      end select

      ! The storage_new may reallocate the p_Rdescriptors array, so get the
      ! pointer again to be sure it is correct and not pointing to nowhere!
      p_rsource => p_rheap%p_Rdescriptors(h_source)
      p_rdest => p_rheap%p_Rdescriptors(h_dest)

    else
      
      ! Check if the given destination handle is compatible with the source handle
      p_rdest => p_rheap%p_Rdescriptors(h_dest)
      
      ! Check if source and destination handle have the same dimension
      if (p_rsource%idimension .ne. p_rdest%idimension) then
        call output_line ('Dimension of source and destination handles are different!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
      end if

    end if

    ! What is to copy
    select case (p_rsource%idataType)
    case (ST_SINGLE)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Fsingle3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Fsingle3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Fsingle3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Fsingle3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Fsingle3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Fsingle3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Fsingle3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Fsingle3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Fsingle3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Fsingle3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_DOUBLE)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Fsingle3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Ddouble3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Ddouble3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Ddouble3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Fsingle3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Ddouble3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Ddouble3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Ddouble3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Ddouble3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Fsingle3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_QUAD)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Fsingle3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Qquad3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Qquad3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Qquad3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Fsingle3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Qquad3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Qquad3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Qquad3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Qquad3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Fsingle3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do

      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Fsingle3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iinteger3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Fsingle3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iinteger3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iinteger3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                p_rsource%p_Fsingle3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Fsingle3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint8_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Fsingle3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint8_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint8_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Fsingle3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Fsingle3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint16_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Fsingle3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint16_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint16_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Fsingle3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Fsingle3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint32_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Fsingle3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint32_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint32_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Fsingle3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Fsingle3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Fsingle3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Fsingle3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint64_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Fsingle3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Fsingle3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Fsingle3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint64_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint64_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Fsingle3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
        
      end select
      
    case (ST_DOUBLE)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Ddouble3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Fsingle3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Fsingle3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Fsingle3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Ddouble3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Fsingle3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Fsingle3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Fsingle3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Fsingle3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Ddouble3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_DOUBLE)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Ddouble3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Ddouble3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Ddouble3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Ddouble3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Ddouble3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Ddouble3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Ddouble3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Ddouble3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Ddouble3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Ddouble3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_QUAD)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Ddouble3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Qquad3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Qquad3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Qquad3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Ddouble3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Qquad3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Qquad3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Qquad3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Qquad3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Ddouble3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Ddouble3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iinteger3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Ddouble3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iinteger3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iinteger3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Ddouble3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Ddouble3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint8_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Ddouble3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint8_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint8_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Ddouble3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Ddouble3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint16_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Ddouble3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint16_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint16_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Ddouble3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Ddouble3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint32_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Ddouble3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint32_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint32_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Ddouble3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do

      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Ddouble3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Ddouble3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Ddouble3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint64_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Ddouble3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Ddouble3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Ddouble3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint64_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint64_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Ddouble3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
        
      end select
      
    case (ST_QUAD)
      select case (p_rdest%idataType)
      case (ST_SINGLE)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Qquad3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Fsingle3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Fsingle3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Fsingle3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Qquad3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Fsingle3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Fsingle3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Fsingle3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Fsingle3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Qquad3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_DOUBLE)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Qquad3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Ddouble3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Ddouble3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Ddouble3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Qquad3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Ddouble3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Ddouble3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Ddouble3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Ddouble3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Qquad3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_QUAD)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Qquad3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Qquad3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Qquad3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Qquad3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Qquad3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Qquad3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Qquad3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Qquad3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Qquad3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Qquad3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Qquad3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iinteger3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Qquad3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iinteger3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iinteger3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Qquad3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Qquad3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint8_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Qquad3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint8_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint8_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Qquad3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Qquad3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint16_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Qquad3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint16_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint16_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Qquad3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Qquad3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint32_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Qquad3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint32_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint32_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Qquad3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do

      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Qquad3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Qquad3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Qquad3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint64_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Qquad3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Qquad3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Qquad3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint64_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint64_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Qquad3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do

      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
        
      end select
      
    case (ST_INT)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iinteger3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iinteger3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iinteger3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iinteger3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iinteger3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iinteger3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iinteger3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint8_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iinteger3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint8_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint8_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iinteger3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iinteger3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint16_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iinteger3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint16_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint16_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iinteger3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iinteger3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint32_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iinteger3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint32_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint32_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iinteger3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iinteger3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iinteger3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iinteger3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint64_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iinteger3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iinteger3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iinteger3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint64_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint64_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iinteger3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
        
      end select
      
    case (ST_INT8)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint8_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iinteger3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint8_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iinteger3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iinteger3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint8_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint8_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint8_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint8_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint8_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint8_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint8_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint8_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint16_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint8_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint16_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint16_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint8_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint8_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint32_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint8_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint32_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint32_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint8_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iint8_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint8_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint8_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint64_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint8_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint8_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint8_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint64_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint64_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint8_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
        
      end select
      
    case (ST_INT16)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint16_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iinteger3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint16_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iinteger3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iinteger3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint16_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint16_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint8_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint16_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint8_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint8_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint16_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint16_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint16_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint16_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint16_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint16_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint16_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint16_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint32_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint16_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint32_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint32_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint16_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iint16_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint16_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint16_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint64_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint16_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint16_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint16_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint64_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint64_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint16_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
        
      end select
      
    case (ST_INT32)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint32_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iinteger3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint32_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iinteger3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iinteger3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint32_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint32_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint8_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint32_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint8_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint8_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint32_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint32_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint16_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint32_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint16_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint16_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint32_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint32_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint32_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint32_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint32_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint32_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint32_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iint32_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint32_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint32_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint64_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint32_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint32_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint32_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint64_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint64_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint32_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
        
      end select
      
    case (ST_INT64)
      select case(p_rdest%idataType)
      case (ST_INT)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint64_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iinteger3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint64_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iinteger3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iinteger3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iinteger3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iinteger3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint64_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT8)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint64_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint8_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint64_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint8_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint8_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint8_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint8_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint64_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT16)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint64_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint16_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint64_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint16_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint16_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint16_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint16_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint64_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT32)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint64_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint32_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint64_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint32_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint32_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint32_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint32_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint64_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case (ST_INT64)
        if (Istart_source(1) < lbound(p_rsource%p_Iint64_3D,1) .or.&
            Istart_source(2) < lbound(p_rsource%p_Iint64_3D,2) .or.&
            Istart_source(3) < lbound(p_rsource%p_Iint64_3D,3) .or.&
            Istart_dest(1)   < lbound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)   < lbound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)   < lbound(p_rdest%p_Iint64_3D,3) .or. &
            Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Iint64_3D,1) .or. &
            Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Iint64_3D,2) .or. &
            Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Iint64_3D,3) .or. &
            Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Iint64_3D,1) .or. &
            Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Iint64_3D,2) .or. &
            Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Iint64_3D,3)) then
          call output_line ('Subarrays incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
          call sys_halt()
        end if
        ! Copy by hand
        do k=1,Ilength(3)
          do j=1,Ilength(2)
            do i=1,Ilength(1)
              p_rdest%p_Iint64_3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                  p_rsource%p_Iint64_3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
            end do
          end do
        end do
        
      case default
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
        
      end select


    case (ST_LOGICAL)
      if (p_rdest%idataType .ne. ST_LOGICAL) then
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
      end if
      if (Istart_source(1) < lbound(p_rsource%p_Blogical3D,1) .or.&
          Istart_source(2) < lbound(p_rsource%p_Blogical3D,2) .or.&
          Istart_source(3) < lbound(p_rsource%p_Blogical3D,3) .or.&
          Istart_dest(1)   < lbound(p_rdest%p_Blogical3D,1) .or. &
          Istart_dest(2)   < lbound(p_rdest%p_Blogical3D,2) .or. &
          Istart_dest(3)   < lbound(p_rdest%p_Blogical3D,3) .or. &
          Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Blogical3D,1) .or. &
          Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Blogical3D,2) .or. &
          Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Blogical3D,3) .or. &
          Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Blogical3D,1) .or. &
          Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Blogical3D,2) .or. &
          Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Blogical3D,3)) then
        call output_line ('Subarrays incompatible!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
      end if
      ! Copy by hand
      do k=1,Ilength(3)
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Blogical3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                p_rsource%p_Blogical3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
          end do
        end do
      end do
      
    case (ST_CHAR)
      if (p_rdest%idataType .ne. ST_CHAR) then
        call output_line ('Unsupported data type combination!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
      end if
      if (Istart_source(1) < lbound(p_rsource%p_Schar3D,1) .or.&
          Istart_source(2) < lbound(p_rsource%p_Schar3D,2) .or.&
          Istart_source(3) < lbound(p_rsource%p_Schar3D,3) .or.&
          Istart_dest(1)   < lbound(p_rdest%p_Schar3D,1) .or. &
          Istart_dest(2)   < lbound(p_rdest%p_Schar3D,2) .or. &
          Istart_dest(3)   < lbound(p_rdest%p_Schar3D,3) .or. &
          Istart_source(1)+Ilength(1)-1 > ubound(p_rsource%p_Schar3D,1) .or. &
          Istart_source(2)+Ilength(2)-1 > ubound(p_rsource%p_Schar3D,2) .or. &
          Istart_source(3)+Ilength(3)-1 > ubound(p_rsource%p_Schar3D,3) .or. &
          Istart_dest(1)+Ilength(1)-1   > ubound(p_rdest%p_Schar3D,1) .or. &
          Istart_dest(2)+Ilength(2)-1   > ubound(p_rdest%p_Schar3D,2) .or. &
          Istart_dest(3)+Ilength(3)-1   > ubound(p_rdest%p_Schar3D,3)) then
        call output_line ('Subarrays incompatible!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
        call sys_halt()
      end if
      ! Copy by hand
      do k=1,Ilength(3)
        do j=1,Ilength(2)
          do i=1,Ilength(1)
            p_rdest%p_Schar3D(Istart_dest(1)+i-1,Istart_dest(2)+j-1,Istart_dest(3)+k-1) = &
                p_rsource%p_Schar3D(Istart_source(1)+i-1,Istart_source(2)+j-1,Istart_source(3)+k-1)
          end do
        end do
      end do
      
    case default
      call output_line ('Unknown data type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_copy_explicit3D')
      call sys_halt()
    end select

  end subroutine storage_copy_explicit3D

!************************************************************************

!<subroutine>

  subroutine storage_info (bprintHandles, rheap)

!<description>
  ! This routine prints information about the current memory consumption
  ! in a memory block to screen.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the handles still remaining in the
  ! heap together with their names are printed to the terminal.
  logical, intent(in), optional :: bprintHandles

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
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
                   trim(sys_smemL(p_rheap%p_Rdescriptors(i)%imemBytes)) //&
                   ', Type=' // trim(sys_siL(p_rheap%p_Rdescriptors(i)%idataType,15)) //&
                   ' Name=' // trim(adjustl(p_rheap%p_Rdescriptors(i)%sname)) )
            else
              call output_line ( &
                   'Handle ' // trim(sys_siL(i,10)) // ', 2D, Length=' // &
                   trim(sys_smemL(p_rheap%p_Rdescriptors(i)%imemBytes)) // &
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
    call output_line ('Memory in use:                   '//&
                      trim(sys_smemL(p_rheap%itotalMem)))
    call output_line ('Current total number of handles: '//&
                      trim(sys_siL(size(p_rheap%p_IfreeHandles),15)))
    call output_line ('Maximum number of handles used:  '//&
                      trim(sys_siL(p_rheap%nhandlesInUseMax,15)))
    call output_line ('Maximum used memory:             '//&
                      trim(sys_smemL(p_rheap%itotalMemMax)))
                      
  end subroutine storage_info

!************************************************************************

!<subroutine>

  subroutine storage_getdatatype (ihandle, idatatype, rheap)

!<description>
  ! Returns the datatype of an array identified by ihandle.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<output>
  ! Datatype of the array identified by ihandle.
  integer, intent(out) :: idatatype
!</output>

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
    case (ST_QUAD)
      idatatype = ST_QUAD
    case (ST_INT)
      idatatype = ST_INT
    case (ST_INT8)
      idatatype = ST_INT8
    case (ST_INT16)
      idatatype = ST_INT16
    case (ST_INT32)
      idatatype = ST_INT32
    case (ST_INT64)
      idatatype = ST_INT64
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
    case default
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
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<output>
  ! Dimension of the array identified by ihandle.
  integer, intent(out) :: idimension
!</output>

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

  subroutine storage_realloc (scall, isize, ihandle, cinitNewBlock, bcopy, rheap)

!<description>
  ! This routine reallocates an existing memory block with a new desired
  ! size. In case of a multiple-dimension block, the last dimension
  ! is changed. isize is the size of the new memory block / the new size
  ! of the last dimension.
  !
  ! Warning: Reallocation of an array destroys all pointers associated with
  ! the corresponding handle!
!</description>

!<input>

  ! Name of the calling routine
  character(LEN=*), intent(in) :: scall

  ! Requested storage size for the memory block / the new size of the last
  ! dimension in the memory block identified by ihandle
  integer, intent(in) :: isize

  ! Init new storage block identifier (ST_NEWBLOCK_ZERO,
  ! ST_NEWBLOCK_NOINIT, ST_NEWBLOCK_ORDERED).
  ! Specifies how to initialise memory if isize > original array size.
  integer, intent(in) :: cinitNewBlock

  ! OPTIONAL: Copy old data.
  ! =TRUE: Copy data of old array to the new one.
  ! =FALSE: Reallocate memory, do not copy old data.
  ! If not specified, TRUE is assumed.
  logical, intent(in), optional :: bcopy

!</input>

!<inputoutput>

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap

  ! Handle of the memory block.
  integer, intent(inout) :: ihandle

!</inputoutput>

!</subroutine>

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode

  ! New storage node
  type(t_storageNode) :: rstorageNode

  ! size and bounds of the old 1-dimensional array
  integer :: isizeOld,ilbound,iubound

  ! size of the old 2-dimensional array
  integer, dimension(2) :: Isize2Dold,Ilbound2D,Iubound2D

  ! size of the old 3-dimensional array
  integer, dimension(3) :: Isize3Dold,Ilbound3D,Iubound3D

  integer :: ier
  logical :: bcopyData

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
      call sys_halt()
    end if

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
        isizeOld = size(p_rnode%p_Fsingle1D)
        ilbound  = lbound(p_rnode%p_Fsingle1D,1)
        iubound  = ubound(p_rnode%p_Fsingle1D,1)
      case (ST_DOUBLE)
        isizeOld = size(p_rnode%p_Ddouble1D)
        ilbound  = lbound(p_rnode%p_Ddouble1D,1)
        iubound  = ubound(p_rnode%p_Ddouble1D,1)
      case (ST_QUAD)
        isizeOld = size(p_rnode%p_Qquad1D)
        ilbound  = lbound(p_rnode%p_Qquad1D,1)
        iubound  = ubound(p_rnode%p_Qquad1D,1)
      case (ST_INT)
        isizeOld = size(p_rnode%p_Iinteger1D)
        ilbound  = lbound(p_rnode%p_Iinteger1D,1)
        iubound  = ubound(p_rnode%p_Iinteger1D,1)
      case (ST_INT8)
        isizeOld = size(p_rnode%p_Iint8_1D)
        ilbound  = lbound(p_rnode%p_Iint8_1D,1)
        iubound  = ubound(p_rnode%p_Iint8_1D,1)
      case (ST_INT16)
        isizeOld = size(p_rnode%p_Iint16_1D)
        ilbound  = lbound(p_rnode%p_Iint16_1D,1)
        iubound  = ubound(p_rnode%p_Iint16_1D,1)
      case (ST_INT32)
        isizeOld = size(p_rnode%p_Iint32_1D)
        ilbound  = lbound(p_rnode%p_Iint32_1D,1)
        iubound  = ubound(p_rnode%p_Iint32_1D,1)
      case (ST_INT64)
        isizeOld = size(p_rnode%p_Iint64_1D)
        ilbound  = lbound(p_rnode%p_Iint64_1D,1)
        iubound  = ubound(p_rnode%p_Iint64_1D,1)
      case (ST_LOGICAL)
        isizeOld = size(p_rnode%p_Blogical1D)
        ilbound  = lbound(p_rnode%p_Blogical1D,1)
        iubound  = ubound(p_rnode%p_Blogical1D,1)
      case (ST_CHAR)
        isizeOld = size(p_rnode%p_Schar1D)
        ilbound  = lbound(p_rnode%p_Schar1D,1)
        iubound  = ubound(p_rnode%p_Schar1D,1)
      end select

      ! Do we really have to change anything?
      if (isize .eq. isizeOld) return

#ifdef USE_C_PTR_STORAGE
      if (ilbound .eq. 1) then
        ! Allocate new memory and initialise it - if it is larger than the old
        ! memory block.
        select case (rstorageNode%idataType)
        case (ST_SINGLE)
          rstorageNode%imemBytes = int(isize,I64)*ST_SINGLE_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Fsingle1D, (/isize/))
        case (ST_DOUBLE)
          rstorageNode%imemBytes = int(isize,I64)*ST_DOUBLE_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Ddouble1D, (/isize/))
        case (ST_QUAD)
          rstorageNode%imemBytes = int(isize,I64)*ST_QUAD_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Qquad1D, (/isize/))
        case (ST_INT)
          rstorageNode%imemBytes = int(isize,I64)*ST_INT_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iinteger1D, (/isize/))
        case (ST_INT8)
          rstorageNode%imemBytes = int(isize,I64)*ST_INT8_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint8_1D, (/isize/))
        case (ST_INT16)
          rstorageNode%imemBytes = int(isize,I64)*ST_INT16_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint16_1D, (/isize/))
        case (ST_INT32)
          rstorageNode%imemBytes = int(isize,I64)*ST_INT32_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint32_1D, (/isize/))
        case (ST_INT64)
          rstorageNode%imemBytes = int(isize,I64)*ST_INT64_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint64_1D, (/isize/))
        case (ST_LOGICAL)
          rstorageNode%imemBytes = int(isize,I64)*ST_LOGICAL_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Blogical1D, (/isize/))
        case (ST_CHAR)
          rstorageNode%imemBytes = int(isize,I64)*ST_CHAR_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Schar1D, (/isize/))
        case default
          call output_line ('Unsupported memory type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
          call sys_halt()
        end select
        
        ! Nasty trick but quick!
        goto 100
      else
        call output_line ('Resorting to standard ALLOCATE for fixed-size memory!', &
                          OU_CLASS_MSG,OU_MODE_STD,'storage_realloc')
      end if
#endif

      ! Allocate new memory and initialise it - if it is larger than the old
      ! memory block.
      select case (rstorageNode%idataType)
      case (ST_SINGLE)
        allocate(rstorageNode%p_Fsingle1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_SINGLE_BYTES
      case (ST_DOUBLE)
        allocate(rstorageNode%p_Ddouble1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_DOUBLE_BYTES
      case (ST_QUAD)
        allocate(rstorageNode%p_Qquad1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_QUAD_BYTES
      case (ST_INT)
        allocate(rstorageNode%p_Iinteger1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_INT_BYTES
      case (ST_INT8)
        allocate(rstorageNode%p_Iint8_1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_INT8_BYTES
      case (ST_INT16)
        allocate(rstorageNode%p_Iint16_1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_INT16_BYTES
      case (ST_INT32)
        allocate(rstorageNode%p_Iint32_1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_INT32_BYTES
      case (ST_INT64)
        allocate(rstorageNode%p_Iint64_1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_INT64_BYTES
      case (ST_LOGICAL)
        allocate(rstorageNode%p_Blogical1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_LOGICAL_BYTES
      case (ST_CHAR)
        allocate(rstorageNode%p_Schar1D(ilbound:ilbound+isize-1))
        rstorageNode%imemBytes = int(isize,I64)*ST_CHAR_BYTES
      case default
        call output_line ('Unsupported memory type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
        call sys_halt()
      end select

100   iubound=min(ilbound+isize-1,iubound)
      if (isize > isizeOld) &
        call storage_initialiseNode (rstorageNode,cinitNewBlock,iubound+1_I32)

      ! Copy old data?
      if (bcopyData) then
        select case (rstorageNode%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorSngl (p_rnode%p_Fsingle1D(ilbound:iubound),&
                                    rstorageNode%p_Fsingle1D(ilbound:iubound))
        case (ST_DOUBLE)
          call lalg_copyVectorDble (p_rnode%p_Ddouble1D(ilbound:iubound),&
                                    rstorageNode%p_Ddouble1D(ilbound:iubound))
        case (ST_QUAD)
          call lalg_copyVectorQuad (p_rnode%p_Qquad1D(ilbound:iubound),&
                                    rstorageNode%p_Qquad1D(ilbound:iubound))
        case (ST_INT)
          call lalg_copyVectorInt (p_rnode%p_Iinteger1D(ilbound:iubound),&
                                   rstorageNode%p_Iinteger1D(ilbound:iubound))
        case (ST_INT8)
          call lalg_copyVectorInt (p_rnode%p_Iint8_1D(ilbound:iubound),&
                                   rstorageNode%p_Iint8_1D(ilbound:iubound))
        case (ST_INT16)
          call lalg_copyVectorInt (p_rnode%p_Iint16_1D(ilbound:iubound),&
                                   rstorageNode%p_Iint16_1D(ilbound:iubound))
        case (ST_INT32)
          call lalg_copyVectorInt (p_rnode%p_Iint32_1D(ilbound:iubound),&
                                   rstorageNode%p_Iint32_1D(ilbound:iubound))
        case (ST_INT64)
          call lalg_copyVectorInt (p_rnode%p_Iint64_1D(ilbound:iubound),&
                                   rstorageNode%p_Iint64_1D(ilbound:iubound))
        case (ST_LOGICAL)
          call lalg_copyVectorLogical (p_rnode%p_Blogical1D(ilbound:iubound),&
                                       rstorageNode%p_Blogical1D(ilbound:iubound))
        case (ST_CHAR)
          call lalg_copyVectorChar (p_rnode%p_Schar1D(ilbound:iubound),&
                                    rstorageNode%p_Schar1D(ilbound:iubound))
        end select
      end if

    case (2)

      ! Get the size of the old storage node.
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        Isize2Dold = shape(p_rnode%p_Fsingle2D)
        Ilbound2D  = lbound(p_rnode%p_Fsingle2D)
        Iubound2D  = ubound(p_rnode%p_Fsingle2D)
      case (ST_DOUBLE)
        Isize2Dold = shape(p_rnode%p_Ddouble2D)
        Ilbound2D  = lbound(p_rnode%p_Ddouble2D)
        Iubound2D  = ubound(p_rnode%p_Ddouble2D)
      case (ST_QUAD)
        Isize2Dold = shape(p_rnode%p_Qquad2D)
        Ilbound2D  = lbound(p_rnode%p_Qquad2D)
        Iubound2D  = ubound(p_rnode%p_Qquad2D)
      case (ST_INT)
        Isize2Dold = shape(p_rnode%p_Iinteger2D)
        Ilbound2D  = lbound(p_rnode%p_Iinteger2D)
        Iubound2D  = ubound(p_rnode%p_Iinteger2D)
      case (ST_INT8)
        Isize2Dold = shape(p_rnode%p_Iint8_2D)
        Ilbound2D  = lbound(p_rnode%p_Iint8_2D)
        Iubound2D  = ubound(p_rnode%p_Iint8_2D)
      case (ST_INT16)
        Isize2Dold = shape(p_rnode%p_Iint16_2D)
        Ilbound2D  = lbound(p_rnode%p_Iint16_2D)
        Iubound2D  = ubound(p_rnode%p_Iint16_2D)
      case (ST_INT32)
        Isize2Dold = shape(p_rnode%p_Iint32_2D)
        Ilbound2D  = lbound(p_rnode%p_Iint32_2D)
        Iubound2D  = ubound(p_rnode%p_Iint32_2D)
      case (ST_INT64)
        Isize2Dold = shape(p_rnode%p_Iint64_2D)
        Ilbound2D  = lbound(p_rnode%p_Iint64_2D)
        Iubound2D  = ubound(p_rnode%p_Iint64_2D)
      case (ST_LOGICAL)
        Isize2Dold = shape(p_rnode%p_Blogical2D)
        Ilbound2D  = lbound(p_rnode%p_Blogical2D)
        Iubound2D  = ubound(p_rnode%p_Blogical2D)
      case (ST_CHAR)
        Isize2Dold = shape(p_rnode%p_Schar2D)
        Ilbound2D  = lbound(p_rnode%p_Schar2D)
        Iubound2D  = ubound(p_rnode%p_Schar2D)
      end select

      ! Do we really have to change anything?
      if (isize .eq. Isize2Dold(2)) return

#ifdef USE_C_PTR_STORAGE
      if (all(Ilbound2D .eq. 1)) then
        ! Allocate new memory and initialise it - if it is larger than the old
        ! memory block.
        select case (rstorageNode%idataType)
        case (ST_SINGLE)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_SINGLE_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Fsingle2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_DOUBLE)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_DOUBLE_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Ddouble2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_QUAD)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_QUAD_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Qquad2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_INT)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iinteger2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_INT8)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT8_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint8_2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_INT16)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT16_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint16_2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_INT32)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT32_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint32_2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_INT64)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT64_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint64_2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_LOGICAL)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_LOGICAL_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Blogical2D,&
                           (/Isize2Dold(1),isize/))
        case (ST_CHAR)
          rstorageNode%imemBytes = int(Isize2Dold(1),I64)*int(isize,I64)*ST_CHAR_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Schar2D,&
                           (/Isize2Dold(1),isize/))
        case default
          call output_line ('Unsupported memory type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
          call sys_halt()
        end select
        
        ! Nasty trick but quick
        goto 200
      else
        call output_line ('Resorting to standard ALLOCATE for fixed-size memory!', &
                          OU_CLASS_MSG,OU_MODE_STD,'storage_realloc')
      end if
#endif

      ! Allocate new memory and initialise it - if it is larger than the old
      ! memory block.
      select case (rstorageNode%idataType)
      case (ST_SINGLE)
        allocate(rstorageNode%p_Fsingle2D(Ilbound2D(1):Iubound2D(1),&
                                          Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_SINGLE_BYTES
      case (ST_DOUBLE)
        allocate(rstorageNode%p_Ddouble2D(Ilbound2D(1):Iubound2D(1),&
                                          Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_DOUBLE_BYTES
      case (ST_QUAD)
        allocate(rstorageNode%p_Qquad2D(Ilbound2D(1):Iubound2D(1),&
                                        Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_QUAD_BYTES
      case (ST_INT)
        allocate(rstorageNode%p_Iinteger2D(Ilbound2D(1):Iubound2D(1),&
                                           Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT_BYTES
      case (ST_INT8)
        allocate(rstorageNode%p_Iint8_2D(Ilbound2D(1):Iubound2D(1),&
                                         Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT8_BYTES
      case (ST_INT16)
        allocate(rstorageNode%p_Iint16_2D(Ilbound2D(1):Iubound2D(1),&
                                          Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT16_BYTES
      case (ST_INT32)
        allocate(rstorageNode%p_Iint32_2D(Ilbound2D(1):Iubound2D(1),&
                                          Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT32_BYTES
      case (ST_INT64)
        allocate(rstorageNode%p_Iint64_2D(Ilbound2D(1):Iubound2D(1),&
                                          Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_INT64_BYTES
      case (ST_LOGICAL)
        allocate(rstorageNode%p_Blogical2D(Ilbound2D(1):Iubound2D(1),&
                                           Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_LOGICAL_BYTES
      case (ST_CHAR)
        allocate(rstorageNode%p_Schar2D(Ilbound2D(1):Iubound2D(1),&
                                        Ilbound2D(2):Ilbound2D(2)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize2Dold(1),I64)*int(isize,I64)*ST_CHAR_BYTES
      case default
        call output_line ('Unsupported memory type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
        call sys_halt()
      end select

200   Iubound2D(2)=min(Ilbound2D(2)+isize-1,Iubound2D(2))
      if (isize > Isize2Dold(2)) &
          call storage_initialiseNode (rstorageNode,cinitNewBlock,&
                                       Iubound2D(2)+1_I32)

      ! Copy old data?
      if (bcopyData) then
        ! Here it is easier than in storage_copy as we can be sure, source and
        ! destination array have the same type!
        select case (rstorageNode%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorSngl2D(p_rnode%p_Fsingle2D(:,:Iubound2D(2)),&
                                     rstorageNode%p_Fsingle2D(:,:Iubound2D(2)))

        case (ST_DOUBLE)
          call lalg_copyVectorDble2D(p_rnode%p_Ddouble2D(:,:Iubound2D(2)),&
                                     rstorageNode%p_Ddouble2D(:,:Iubound2D(2)))

        case (ST_QUAD)
          call lalg_copyVectorQuad2D(p_rnode%p_Qquad2D(:,:Iubound2D(2)),&
                                     rstorageNode%p_Qquad2D(:,:Iubound2D(2)))

        case (ST_INT)
          call lalg_copyVectorInt2D(p_rnode%p_Iinteger2D(:,:Iubound2D(2)),&
                                    rstorageNode%p_Iinteger2D(:,:Iubound2D(2)))

        case (ST_INT8)
          call lalg_copyVectorInt2D(p_rnode%p_Iint8_2D(:,:Iubound2D(2)),&
                                    rstorageNode%p_Iint8_2D(:,:Iubound2D(2)))

        case (ST_INT16)
          call lalg_copyVectorInt2D(p_rnode%p_Iint16_2D(:,:Iubound2D(2)),&
                                    rstorageNode%p_Iint16_2D(:,:Iubound2D(2)))

        case (ST_INT32)
          call lalg_copyVectorInt2D(p_rnode%p_Iint32_2D(:,:Iubound2D(2)),&
                                    rstorageNode%p_Iint32_2D(:,:Iubound2D(2)))

        case (ST_INT64)
          call lalg_copyVectorInt2D(p_rnode%p_Iint64_2D(:,:Iubound2D(2)),&
                                    rstorageNode%p_Iint64_2D(:,:Iubound2D(2)))

        case (ST_LOGICAL)
          call lalg_copyVectorLogical2D(p_rnode%p_Blogical2D(:,:Iubound2D(2)),&
                                        rstorageNode%p_Blogical2D(:,:Iubound2D(2)))

        case (ST_CHAR)
          call lalg_copyVectorChar2D(p_rnode%p_Schar2D(:,:Iubound2D(2)),&
                                     rstorageNode%p_Schar2D(:,:Iubound2D(2)))
        end select

      end if

    case (3)

      ! Get the size of the old storage node.
      select case (p_rnode%idataType)
      case (ST_SINGLE)
        Isize3Dold = shape(p_rnode%p_Fsingle3D)
        ilbound3D  = lbound(p_rnode%p_Fsingle3D)
        Iubound3D  = ubound(p_rnode%p_Fsingle3D)
      case (ST_DOUBLE)
        Isize3Dold = shape(p_rnode%p_Ddouble3D)
        ilbound3D  = lbound(p_rnode%p_Ddouble3D)
        Iubound3D  = ubound(p_rnode%p_Ddouble3D)
      case (ST_QUAD)
        Isize3Dold = shape(p_rnode%p_Qquad3D)
        ilbound3D  = lbound(p_rnode%p_Qquad3D)
        Iubound3D  = ubound(p_rnode%p_Qquad3D)
      case (ST_INT)
        Isize3Dold = shape(p_rnode%p_Iinteger3D)
        ilbound3D  = lbound(p_rnode%p_Iinteger3D)
        Iubound3D  = ubound(p_rnode%p_Iinteger3D)
      case (ST_INT8)
        Isize3Dold = shape(p_rnode%p_Iint8_3D)
        ilbound3D  = lbound(p_rnode%p_Iint8_3D)
        Iubound3D  = ubound(p_rnode%p_Iint8_3D)
      case (ST_INT16)
        Isize3Dold = shape(p_rnode%p_Iint16_3D)
        ilbound3D  = lbound(p_rnode%p_Iint16_3D)
        Iubound3D  = ubound(p_rnode%p_Iint16_3D)
      case (ST_INT32)
        Isize3Dold = shape(p_rnode%p_Iint32_3D)
        ilbound3D  = lbound(p_rnode%p_Iint32_3D)
        Iubound3D  = ubound(p_rnode%p_Iint32_3D)
      case (ST_INT64)
        Isize3Dold = shape(p_rnode%p_Iint64_3D)
        ilbound3D  = lbound(p_rnode%p_Iint64_3D)
        Iubound3D  = ubound(p_rnode%p_Iint64_3D)
      case (ST_LOGICAL)
        Isize3Dold = shape(p_rnode%p_Blogical3D)
        ilbound3D  = lbound(p_rnode%p_Blogical3D)
        Iubound3D  = ubound(p_rnode%p_Blogical3D)
      case (ST_CHAR)
        Isize3Dold = shape(p_rnode%p_Schar3D)
        ilbound3D  = lbound(p_rnode%p_Schar3D)
        Iubound3D  = ubound(p_rnode%p_Schar3D)
      end select

      ! Do we really have to change anything?
      if (isize .eq. Isize3Dold(3)) return

#ifdef USE_C_PTR_STORAGE
      if (all(Ilbound3D .eq. 1)) then
        ! Allocate new memory and initialise it - if it is larger than the old
        ! memory block.
        select case (rstorageNode%idataType)
        case (ST_SINGLE)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_SINGLE_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Fsingle3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_DOUBLE)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_DOUBLE_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Ddouble3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_QUAD)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_QUAD_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Qquad3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_INT)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_INT_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iinteger3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_INT8)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_INT8_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint8_3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_INT16)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_INT16_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint16_3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_INT32)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_INT32_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint32_3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_INT64)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_INT64_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Iint64_3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_LOGICAL)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_LOGICAL_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Blogical3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case (ST_CHAR)
          rstorageNode%imemBytes = int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
                                   int(isize,I64)*ST_CHAR_BYTES
          ier = c_allocate(rstorageNode%chostMemPtr, rstorageNode%imemBytes)
          call c_f_pointer(rstorageNode%chostMemPtr, rstorageNode%p_Schar3D,&
                           (/Isize3Dold(1),Isize3Dold(2),isize/))
        case default
          call output_line ('Unsupported memory type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
          call sys_halt()
        end select
      
        ! Nasty trick but quick
        goto 300
      else
        call output_line ('Resorting to standard ALLOCATE for fixed-size memory!', &
                          OU_CLASS_MSG,OU_MODE_STD,'storage_realloc')
      end if
#endif
      
      ! Allocate new memory and initialise it - if it is larger than the old
      ! memory block.
      select case (rstorageNode%idataType)
      case (ST_SINGLE)
        allocate(rstorageNode%p_Fsingle3D(Ilbound3D(1):Iubound3D(1),&
                                          Ilbound3D(2):Iubound3D(2),&
                                          Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_SINGLE_BYTES
      case (ST_DOUBLE)
        allocate(rstorageNode%p_Ddouble3D(Ilbound3D(1):Iubound3D(1),&
                                          Ilbound3D(2):Iubound3D(2),&
                                          Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_DOUBLE_BYTES
      case (ST_QUAD)
        allocate(rstorageNode%p_Qquad3D(Ilbound3D(1):Iubound3D(1),&
                                        Ilbound3D(2):Iubound3D(2),&
                                        Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_QUAD_BYTES
      case (ST_INT)
        allocate(rstorageNode%p_Iinteger3D(Ilbound3D(1):Iubound3D(1),&
                                           Ilbound3D(2):Iubound3D(2),&
                                           Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_INT_BYTES
      case (ST_INT8)
        allocate(rstorageNode%p_Iint8_3D(Ilbound3D(1):Iubound3D(1),&
                                         Ilbound3D(2):Iubound3D(2),&
                                         Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_INT8_BYTES
      case (ST_INT16)
        allocate(rstorageNode%p_Iint16_3D(Ilbound3D(1):Iubound3D(1),&
                                          Ilbound3D(2):Iubound3D(2),&
                                          Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_INT16_BYTES
      case (ST_INT32)
        allocate(rstorageNode%p_Iint32_3D(Ilbound3D(1):Iubound3D(1),&
                                          Ilbound3D(2):Iubound3D(2),&
                                          Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_INT32_BYTES
      case (ST_INT64)
        allocate(rstorageNode%p_Iint64_3D(Ilbound3D(1):Iubound3D(1),&
                                          Ilbound3D(2):Iubound3D(2),&
                                          Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_INT64_BYTES
      case (ST_LOGICAL)
        allocate(rstorageNode%p_Blogical3D(Ilbound3D(1):Iubound3D(1),&
                                           Ilbound3D(2):Iubound3D(2),&
                                           Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_LOGICAL_BYTES
      case (ST_CHAR)
        allocate(rstorageNode%p_Schar3D(Ilbound3D(1):Iubound3D(1),&
                                        Ilbound3D(2):Iubound3D(2),&
                                        Ilbound3D(3):Ilbound3D(3)+isize-1))
        rstorageNode%imemBytes = &
             int(Isize3Dold(1),I64)*int(Isize3Dold(2),I64)*&
             int(isize,I64)*ST_CHAR_BYTES
      case default
        call output_line ('Unsupported memory type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
        call sys_halt()
      end select

300   Iubound3D(3)=min(Ilbound3D(3)+isize-1,Iubound3D(3))
      if (isize > Isize3Dold(3)) &
          call storage_initialiseNode (rstorageNode,cinitNewBlock,&
                                       Iubound3D(3)+1_I32)

      ! Copy old data?
      if (bcopyData) then
        ! Here it is easier than in storage_copy as we can be sure, source and
        ! destination array have the same type!
        select case (rstorageNode%idataType)
        case (ST_SINGLE)
          call lalg_copyVectorSngl3D(p_rnode%p_Fsingle3D(:,:,:Iubound3D(3)),&
                                     rstorageNode%p_Fsingle3D(:,:,:Iubound3D(3)))

        case (ST_DOUBLE)
          call lalg_copyVectorDble3D(p_rnode%p_Ddouble3D(:,:,:Iubound3D(3)),&
                                     rstorageNode%p_Ddouble3D(:,:,:Iubound3D(3)))

        case (ST_QUAD)
          call lalg_copyVectorQuad3D(p_rnode%p_Qquad3D(:,:,:Iubound3D(3)),&
                                     rstorageNode%p_Qquad3D(:,:,:Iubound3D(3)))

        case (ST_INT)
          call lalg_copyVectorInt3D(p_rnode%p_Iinteger3D(:,:,:Iubound3D(3)),&
                                    rstorageNode%p_Iinteger3D(:,:,:Iubound3D(3)))

        case (ST_INT8)
          call lalg_copyVectorInt3D(p_rnode%p_Iint8_3D(:,:,:Iubound3D(3)),&
                                    rstorageNode%p_Iint8_3D(:,:,:Iubound3D(3)))

        case (ST_INT16)
          call lalg_copyVectorInt3D(p_rnode%p_Iint16_3D(:,:,:Iubound3D(3)),&
                                    rstorageNode%p_Iint16_3D(:,:,:Iubound3D(3)))

        case (ST_INT32)
          call lalg_copyVectorInt3D(p_rnode%p_Iint32_3D(:,:,:Iubound3D(3)),&
                                    rstorageNode%p_Iint32_3D(:,:,:Iubound3D(3)))

        case (ST_INT64)
          call lalg_copyVectorInt3D(p_rnode%p_Iint64_3D(:,:,:Iubound3D(3)),&
                                    rstorageNode%p_Iint64_3D(:,:,:Iubound3D(3)))

        case (ST_LOGICAL)
          call lalg_copyVectorLogical3D(p_rnode%p_Blogical3D(:,:,:Iubound3D(3)),&
                                        rstorageNode%p_Blogical3D(:,:,:Iubound3D(3)))

        case (ST_CHAR)
          call lalg_copyVectorChar3D(p_rnode%p_Schar3D(:,:,:Iubound3D(3)),&
                                     rstorageNode%p_Schar3D(:,:,:Iubound3D(3)))
        end select

      end if

    case default
      call output_line ('Handle '//trim(sys_siL(ihandle,11))//' is neither 1-, 2- nor 3-dimensional!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
      call sys_halt()

    end select

    ! Respect also the temporary memory in the total amount of memory used.
    !$omp critical(storage_global_heap_modify)
    if ((p_rheap%itotalMem + rstorageNode%imemBytes) .gt. p_rheap%itotalMemMax) &
      p_rheap%itotalMemMax = p_rheap%itotalMem + rstorageNode%imemBytes
    !$omp end critical(storage_global_heap_modify)

#ifdef USE_C_PTR_STORAGE
    ! Release host memory physically
    if (storage_isAssociated(p_rnode%chostMemPtr)) then
      if (c_deallocate(p_rnode%chostMemPtr) .gt. 0) then
        call output_line ('Error in freeing memory with c_deallocate!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_realloc')
        call sys_halt()
      end if
    end if

    ! Nullify pointers
    nullify(p_rnode%p_Fsingle1D)
    nullify(p_rnode%p_Ddouble1D)
    nullify(p_rnode%p_Qquad1D)
    nullify(p_rnode%p_Iinteger1D)
    nullify(p_rnode%p_Iint8_1D)
    nullify(p_rnode%p_Iint16_1D)
    nullify(p_rnode%p_Iint32_1D)
    nullify(p_rnode%p_Iint64_1D)
    nullify(p_rnode%p_Blogical1D)
    nullify(p_rnode%p_Schar1D)

    nullify(p_rnode%p_Fsingle2D)
    nullify(p_rnode%p_Ddouble2D)
    nullify(p_rnode%p_Qquad2D)
    nullify(p_rnode%p_Iinteger2D)
    nullify(p_rnode%p_Iint8_2D)
    nullify(p_rnode%p_Iint16_2D)
    nullify(p_rnode%p_Iint32_2D)
    nullify(p_rnode%p_Iint64_2D)
    nullify(p_rnode%p_Blogical2D)
    nullify(p_rnode%p_Schar2D)

    nullify(p_rnode%p_Fsingle3D)
    nullify(p_rnode%p_Ddouble3D)
    nullify(p_rnode%p_Qquad3D)
    nullify(p_rnode%p_Iinteger3D)
    nullify(p_rnode%p_Iint8_3D)
    nullify(p_rnode%p_Iint16_3D)
    nullify(p_rnode%p_Iint32_3D)
    nullify(p_rnode%p_Iint64_3D)
    nullify(p_rnode%p_Blogical3D)
    nullify(p_rnode%p_Schar3D)
#endif

#ifdef ENABLE_COPROCESSOR_SUPPORT
    if (storage_isAssociated(p_rnode%cdeviceMemPtr))&
        call coproc_freeMemoryOnDevice(p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
#endif

    ! Release old data
    if (associated(p_rnode%p_Fsingle1D))  deallocate(p_rnode%p_Fsingle1D)
    if (associated(p_rnode%p_Ddouble1D))  deallocate(p_rnode%p_Ddouble1D)
    if (associated(p_rnode%p_Qquad1D))    deallocate(p_rnode%p_Qquad1D)
    if (associated(p_rnode%p_Iinteger1D)) deallocate(p_rnode%p_Iinteger1D)
    if (associated(p_rnode%p_Iint8_1D))   deallocate(p_rnode%p_Iint8_1D)
    if (associated(p_rnode%p_Iint16_1D))  deallocate(p_rnode%p_Iint16_1D)
    if (associated(p_rnode%p_Iint32_1D))  deallocate(p_rnode%p_Iint32_1D)
    if (associated(p_rnode%p_Iint64_1D))  deallocate(p_rnode%p_Iint64_1D)
    if (associated(p_rnode%p_Blogical1D)) deallocate(p_rnode%p_Blogical1D)
    if (associated(p_rnode%p_Schar1D))    deallocate(p_rnode%p_Schar1D)

    if (associated(p_rnode%p_Fsingle2D))  deallocate(p_rnode%p_Fsingle2D)
    if (associated(p_rnode%p_Ddouble2D))  deallocate(p_rnode%p_Ddouble2D)
    if (associated(p_rnode%p_Qquad2D))    deallocate(p_rnode%p_Qquad2D)
    if (associated(p_rnode%p_Iinteger2D)) deallocate(p_rnode%p_Iinteger2D)
    if (associated(p_rnode%p_Iint8_2D))   deallocate(p_rnode%p_Iint8_2D)
    if (associated(p_rnode%p_Iint16_2D))  deallocate(p_rnode%p_Iint16_2D)
    if (associated(p_rnode%p_Iint32_2D))  deallocate(p_rnode%p_Iint32_2D)
    if (associated(p_rnode%p_Iint64_2D))  deallocate(p_rnode%p_Iint64_2D)
    if (associated(p_rnode%p_Blogical2D)) deallocate(p_rnode%p_Blogical2D)
    if (associated(p_rnode%p_Schar2D))    deallocate(p_rnode%p_Schar2D)

    if (associated(p_rnode%p_Fsingle3D))  deallocate(p_rnode%p_Fsingle3D)
    if (associated(p_rnode%p_Ddouble3D))  deallocate(p_rnode%p_Ddouble3D)
    if (associated(p_rnode%p_Qquad3D))    deallocate(p_rnode%p_Qquad3D)
    if (associated(p_rnode%p_Iinteger3D)) deallocate(p_rnode%p_Iinteger3D)
    if (associated(p_rnode%p_Iint8_3D))   deallocate(p_rnode%p_Iint8_3D)
    if (associated(p_rnode%p_Iint16_3D))  deallocate(p_rnode%p_Iint16_3D)
    if (associated(p_rnode%p_Iint32_3D))  deallocate(p_rnode%p_Iint32_3D)
    if (associated(p_rnode%p_Iint64_3D))  deallocate(p_rnode%p_Iint64_3D)
    if (associated(p_rnode%p_Blogical3D)) deallocate(p_rnode%p_Blogical3D)
    if (associated(p_rnode%p_Schar3D))    deallocate(p_rnode%p_Schar3D)

    ! Correct the memory statistics
    !$omp critical(storage_global_heap_modify)
    p_rheap%itotalMem = p_rheap%itotalMem &
                      - p_rnode%imemBytes + rstorageNode%imemBytes
    if (p_rheap%itotalMem .gt. p_rheap%itotalMemMax) &
      p_rheap%itotalMemMax = p_rheap%itotalMem
    !$omp end critical(storage_global_heap_modify)

    ! Replace the old node by the new one, finish
    p_rnode = rstorageNode

  end subroutine storage_realloc

!************************************************************************

!<function>

  function storage_isEqual (ihandle1, ihandle2, rheap1, rheap2) result (bisequal)

!<description>

  ! This function checks if the content of two different handles is equal

!</description>

!<input>

  ! The first handle
  integer, intent(in) :: ihandle1

  ! The second handle
  integer, intent(in) :: ihandle2

  ! OPTIONAL: first local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap1

  ! OPTIONAL: second local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap2

!</input>

!<result>

  ! .TRUE. if the content is equal.
  logical :: bisequal

!</result>

!</function>

  ! Pointer to the heaps
  type(t_storageBlock), pointer :: p_rheap1,p_rheap2

  ! Identifier for data type
  integer :: idatatype1, idatatype2

  ! Identifier for data dimension
  integer :: idimension1, idimension2

  ! Identifier for size
  integer :: isize1, isize2
  integer, dimension(2) :: Isize2D1, Isize2D2
  integer, dimension(3) :: Isize3D1, Isize3D2

  ! Auxiliary arrays
  real(SP), dimension(:,:,:), pointer     :: p_Fsingle3D1,p_Fsingle3D2
  real(SP), dimension(:,:), pointer       :: p_Fsingle2D1,p_Fsingle2D2
  real(SP), dimension(:),   pointer       :: p_Fsingle1D1,p_Fsingle1D2
  real(DP), dimension(:,:,:), pointer     :: p_Ddouble3D1,p_Ddouble3D2
  real(DP), dimension(:,:), pointer       :: p_Ddouble2D1,p_Ddouble2D2
  real(DP), dimension(:),   pointer       :: p_Ddouble1D1,p_Ddouble1D2
  real(QP), dimension(:,:,:), pointer     :: p_Qquad3D1,p_Qquad3D2
  real(QP), dimension(:,:), pointer       :: p_Qquad2D1,p_Qquad2D2
  real(QP), dimension(:),   pointer       :: p_Qquad1D1,p_Qquad1D2
  integer, dimension(:,:,:), pointer      :: p_Iinteger3D1,p_Iinteger3D2
  integer, dimension(:,:), pointer        :: p_Iinteger2D1,p_Iinteger2D2
  integer, dimension(:),   pointer        :: p_Iinteger1D1,p_Iinteger1D2
  integer(I8), dimension(:,:,:), pointer  :: p_Iint8_3D1,p_Iint8_3D2
  integer(I8), dimension(:,:), pointer    :: p_Iint8_2D1,p_Iint8_2D2
  integer(I8), dimension(:),   pointer    :: p_Iint8_1D1,p_Iint8_1D2
  integer(I16), dimension(:,:,:), pointer :: p_Iint16_3D1,p_Iint16_3D2
  integer(I16), dimension(:,:), pointer   :: p_Iint16_2D1,p_Iint16_2D2
  integer(I16), dimension(:),   pointer   :: p_Iint16_1D1,p_Iint16_1D2
  integer(I32), dimension(:,:,:), pointer :: p_Iint32_3D1,p_Iint32_3D2
  integer(I32), dimension(:,:), pointer   :: p_Iint32_2D1,p_Iint32_2D2
  integer(I32), dimension(:),   pointer   :: p_Iint32_1D1,p_Iint32_1D2
  integer(I64), dimension(:,:,:), pointer :: p_Iint64_3D1,p_Iint64_3D2
  integer(I64), dimension(:,:), pointer   :: p_Iint64_2D1,p_Iint64_2D2
  integer(I64), dimension(:),   pointer   :: p_Iint64_1D1,p_Iint64_1D2
  logical, dimension(:,:,:), pointer      :: p_Blogical3D1,p_Blogical3D2
  logical, dimension(:,:), pointer        :: p_Blogical2D1,p_Blogical2D2
  logical, dimension(:),   pointer        :: p_Blogical1D1,p_Blogical1D2
  character, dimension(:,:,:), pointer    :: p_Schar3D1,p_Schar3D2
  character, dimension(:,:), pointer      :: p_Schar2D1,p_Schar2D2
  character, dimension(:),   pointer      :: p_Schar1D1,p_Schar1D2
  
    ! Get the heaps to use - local or global ones.
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
      case (ST_SINGLE)
        call storage_getbase_single(ihandle1, p_Fsingle1D1, p_rheap1)
        call storage_getbase_single(ihandle2, p_Fsingle1D2, p_rheap2)
        
        bisequal = all(p_Fsingle1D1 .eq. p_Fsingle1D2)

      case (ST_DOUBLE)
        call storage_getbase_double(ihandle1, p_Ddouble1D1, p_rheap1)
        call storage_getbase_double(ihandle2, p_Ddouble1D2, p_rheap2)

        bisequal = all(p_Ddouble1D1 .eq. p_Ddouble1D2)

      case (ST_QUAD)
        call storage_getbase_quad(ihandle1, p_Qquad1D1, p_rheap1)
        call storage_getbase_quad(ihandle2, p_Qquad1D2, p_rheap2)

        bisequal = all(p_Qquad1D1 .eq. p_Qquad1D2)
        
      case (ST_INT)
        call storage_getbase_int(ihandle1, p_Iinteger1D1, p_rheap1)
        call storage_getbase_int(ihandle2, p_Iinteger1D2, p_rheap2)

        bisequal = all(p_Iinteger1D1 .eq. p_Iinteger1D2)
        
      case (ST_INT8)
        call storage_getbase_int8(ihandle1, p_Iint8_1D1, p_rheap1)
        call storage_getbase_int8(ihandle2, p_Iint8_1D2, p_rheap2)

        bisequal = all(p_Iint8_1D1 .eq. p_Iint8_1D2)
        
      case (ST_INT16)
        call storage_getbase_int16(ihandle1, p_Iint16_1D1, p_rheap1)
        call storage_getbase_int16(ihandle2, p_Iint16_1D2, p_rheap2)

        bisequal = all(p_Iint16_1D1 .eq. p_Iint16_1D2)
        
      case (ST_INT32)
        call storage_getbase_int32(ihandle1, p_Iint32_1D1, p_rheap1)
        call storage_getbase_int32(ihandle2, p_Iint32_1D2, p_rheap2)

        bisequal = all(p_Iint32_1D1 .eq. p_Iint32_1D2)
        
      case (ST_INT64)
        call storage_getbase_int64(ihandle1, p_Iint64_1D1, p_rheap1)
        call storage_getbase_int64(ihandle2, p_Iint64_1D2, p_rheap2)

        bisequal = all(p_Iint64_1D1 .eq. p_Iint64_1D2)

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
      call storage_getsize(ihandle1, Isize2D1, p_rheap1)
      call storage_getsize(ihandle2, Isize2D2, p_rheap2)

      if (any(Isize2D1 .ne. Isize2D2)) then
        bisequal = .false.
        return
      end if

      ! What data type do we have?
      select case(idatatype1)
      case (ST_SINGLE)
        call storage_getbase_single2d(ihandle1, p_Fsingle2D1, p_rheap1)
        call storage_getbase_single2d(ihandle2, p_Fsingle2D2, p_rheap2)
        
        bisequal = all(p_Fsingle2D1 .eq. p_Fsingle2D2)

      case (ST_DOUBLE)
        call storage_getbase_double2d(ihandle1, p_Ddouble2D1, p_rheap1)
        call storage_getbase_double2d(ihandle2, p_Ddouble2D2, p_rheap2)

        bisequal = all(p_Ddouble2D1 .eq. p_Ddouble2D2)

      case (ST_QUAD)
        call storage_getbase_quad2d(ihandle1, p_Qquad2D1, p_rheap1)
        call storage_getbase_quad2d(ihandle2, p_Qquad2D2, p_rheap2)

        bisequal = all(p_Qquad2D1 .eq. p_Qquad2D2)

      case (ST_INT)
        call storage_getbase_int2d(ihandle1, p_Iinteger2D1, p_rheap1)
        call storage_getbase_int2d(ihandle2, p_Iinteger2D2, p_rheap2)

        bisequal = all(p_Iinteger2D1 .eq. p_Iinteger2D2)

      case (ST_INT8)
        call storage_getbase_int8_2d(ihandle1, p_Iint8_2D1, p_rheap1)
        call storage_getbase_int8_2d(ihandle2, p_Iint8_2D2, p_rheap2)

        bisequal = all(p_Iint8_2D1 .eq. p_Iint8_2D2)

      case (ST_INT16)
        call storage_getbase_int16_2d(ihandle1, p_Iint16_2D1, p_rheap1)
        call storage_getbase_int16_2d(ihandle2, p_Iint16_2D2, p_rheap2)

        bisequal = all(p_Iint16_2D1 .eq. p_Iint16_2D2)

      case (ST_INT32)
        call storage_getbase_int32_2d(ihandle1, p_Iint32_2D1, p_rheap1)
        call storage_getbase_int32_2d(ihandle2, p_Iint32_2D2, p_rheap2)

        bisequal = all(p_Iint32_2D1 .eq. p_Iint32_2D2)

      case (ST_INT64)
        call storage_getbase_int64_2d(ihandle1, p_Iint64_2D1, p_rheap1)
        call storage_getbase_int64_2d(ihandle2, p_Iint64_2D2, p_rheap2)

        bisequal = all(p_Iint64_2D1 .eq. p_Iint64_2D2)
        
      case (ST_LOGICAL)
        call storage_getbase_logical2d(ihandle1, p_Blogical2D1, p_rheap1)
        call storage_getbase_logical2d(ihandle2, p_Blogical2D2, p_rheap2)

        bisequal = all(p_Blogical2D1 .and. p_Blogical2D2)
        
      case (ST_CHAR)
        call storage_getbase_char2d(ihandle1, p_Schar2D1, p_rheap1)
        call storage_getbase_char2d(ihandle2, p_Schar2D2, p_rheap2)

        bisequal = all(p_Schar2D1 .eq. p_Schar2D2)
      end select

    case (3)

      ! Determine size
      call storage_getsize(ihandle1, Isize3D1, p_rheap1)
      call storage_getsize(ihandle2, Isize3D2, p_rheap2)

      if (any(Isize3D1 .ne. Isize3D2)) then
        bisequal = .false.
        return
      end if

      ! What data type do we have?
      select case(idatatype1)
      case (ST_SINGLE)
        call storage_getbase_single3d(ihandle1, p_Fsingle3D1, p_rheap1)
        call storage_getbase_single3d(ihandle2, p_Fsingle3D2, p_rheap2)
        
        bisequal = all(p_Fsingle3D1 .eq. p_Fsingle3D2)

      case (ST_DOUBLE)
        call storage_getbase_double3d(ihandle1, p_Ddouble3D1, p_rheap1)
        call storage_getbase_double3d(ihandle2, p_Ddouble3D2, p_rheap2)

        bisequal = all(p_Ddouble3D1 .eq. p_Ddouble3D2)

      case (ST_QUAD)
        call storage_getbase_quad3d(ihandle1, p_Qquad3D1, p_rheap1)
        call storage_getbase_quad3d(ihandle2, p_Qquad3D2, p_rheap2)

        bisequal = all(p_Qquad3D1 .eq. p_Qquad3D2)

      case (ST_INT)
        call storage_getbase_int3d(ihandle1, p_Iinteger3D1, p_rheap1)
        call storage_getbase_int3d(ihandle2, p_Iinteger3D2, p_rheap2)

        bisequal = all(p_Iinteger3D1 .eq. p_Iinteger3D2)

      case (ST_INT8)
        call storage_getbase_int8_3d(ihandle1, p_Iint8_3D1, p_rheap1)
        call storage_getbase_int8_3d(ihandle2, p_Iint8_3D2, p_rheap2)

        bisequal = all(p_Iint8_3D1 .eq. p_Iint8_3D2)

      case (ST_INT16)
        call storage_getbase_int16_3d(ihandle1, p_Iint16_3D1, p_rheap1)
        call storage_getbase_int16_3d(ihandle2, p_Iint16_3D2, p_rheap2)

        bisequal = all(p_Iint16_3D1 .eq. p_Iint16_3D2)

      case (ST_INT32)
        call storage_getbase_int32_3d(ihandle1, p_Iint32_3D1, p_rheap1)
        call storage_getbase_int32_3d(ihandle2, p_Iint32_3D2, p_rheap2)

        bisequal = all(p_Iint32_3D1 .eq. p_Iint32_3D2)

      case (ST_INT64)
        call storage_getbase_int64_3d(ihandle1, p_Iint64_3D1, p_rheap1)
        call storage_getbase_int64_3d(ihandle2, p_Iint64_3D2, p_rheap2)

        bisequal = all(p_Iint64_3D1 .eq. p_Iint64_3D2)
        
      case (ST_LOGICAL)
        call storage_getbase_logical3d(ihandle1, p_Blogical3D1, p_rheap1)
        call storage_getbase_logical3d(ihandle2, p_Blogical3D2, p_rheap2)

        bisequal = all(p_Blogical3D1 .and. p_Blogical3D2)
        
      case (ST_CHAR)
        call storage_getbase_char3d(ihandle1, p_Schar3D1, p_rheap1)
        call storage_getbase_char3d(ihandle2, p_Schar3D2, p_rheap2)

        bisequal = all(p_Schar3D1 .eq. p_Schar3D2)
      end select

    case default
      call output_line ('Handle '//trim(sys_siL(ihandle1,11))//' is neither 1-, 2- nor 3-dimensional!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_isEqual')
      call sys_halt()
    end select
    
  end function storage_isEqual

!************************************************************************

!<subroutine>

  subroutine storage_createFpdbObject (rfpdbObjectItem, sname, rheap)

!<description>
  ! This subroutine creates an abstract ObjectItem from the heap
  ! structure that can be stored in the persistence database
!</description>

!<input>
  ! The full qualified name of the object
  character(LEN=*), intent(in) :: sname
  
  ! OPTIONAL: local heap structure to initialise.
  ! If not given, the global heap is initialised.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<inputoutput>
  ! The ObjectItem that is created
  type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,ihandlesInUse

  ! Pointer to the heap to initialise
  type(t_storageBlock), pointer :: p_rheap
  type(t_fpdbDataItem), pointer :: p_fpdbDataItem

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    ! Check if storage block has UUID; otherwise create one
    if (uuid_isNil(p_rheap%ruuid)) then
      call uuid_createUUID(4, p_rheap%ruuid)
    end if
    rfpdbObjectItem%ruuid = p_rheap%ruuid

    ! Set the name and type of the object
    rfpdbObjectItem%sname = sname
    rfpdbObjectItem%stype = 't_storageBlock'
    
    ! Allocate the array of data items: 9 + the number of handles in use
    allocate(rfpdbObjectItem%p_RfpdbDataItem(9+p_rheap%ihandlesInUse))

    ! Fill the array of data items: p_inextFreeHandle
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'p_inextFreeHandle'
    p_fpdbDataItem%iinteger = p_rheap%p_inextFreeHandle

    ! Fill the array of data items: p_ilastFreeHandle
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'p_ilastFreeHandle'
    p_fpdbDataItem%iinteger = p_rheap%p_ilastFreeHandle

    ! Fill the array of data items: ihandlesInUse
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'ihandlesInUse'
    p_fpdbDataItem%iinteger = p_rheap%ihandlesInUse

    ! Fill the array of data items: nhandlesTotal
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'nhandlesTotal'
    p_fpdbDataItem%iinteger = p_rheap%nhandlesTotal

    ! Fill the array of data items: ihandlesDelta
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'ihandlesDelta'
    p_fpdbDataItem%iinteger = p_rheap%ihandlesDelta

    ! Fill the array of data items: itotalMem
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    p_fpdbDataItem%ctype    = FPDB_INT64
    p_fpdbDataItem%sname    = 'itotalMem'
    p_fpdbDataItem%iint64   = p_rheap%itotalMem

    ! Fill the array of data items: nhandlesInUseMax
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'nhandlesInUseMax'
    p_fpdbDataItem%iinteger = p_rheap%nhandlesInUseMax

    ! Fill the array of data items: itotalMemMax
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    p_fpdbDataItem%ctype    = FPDB_INT64
    p_fpdbDataItem%sname    = 'itotalMemMax'
    p_fpdbDataItem%iint64   = p_rheap%itotalMemMax

      ! Fill the array of data items: p_IfreeHandles
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    p_fpdbDataItem%ctype        =  FPDB_INT1D
    p_fpdbDataItem%sname        =  'p_IfreeHandles'
    p_fpdbDataItem%p_Iinteger1D => p_rheap%p_IfreeHandles

    ! Attach all handles currently in use
    ihandlesInUse = 0
    do i = 1, size(p_rheap%p_Rdescriptors)
      
      ! Only attach handles which are in use
      if (p_rheap%p_Rdescriptors(i)%idataType .ne. ST_NOHANDLE) then
        ihandlesInUse = ihandlesInUse+1

        p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9+ihandlesInUse)
        p_fpdbDataItem%ctype = FPDB_OBJECT
        p_fpdbDataItem%sname = 'handle_'//trim(sys_siL(i,11))

        allocate(p_fpdbDataItem%p_fpdbObjectItem)
        call createFpdbObjectHandle(p_fpdbDataItem%p_fpdbObjectItem, i)
      end if
    end do
      
    ! Check consistency of handles in use
    if (ihandlesInUse .ne. p_rheap%ihandlesInUse) then
      call output_line ('Inconsistency detected', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_createFpdbObject')
      call sys_halt()
    end if


  contains

    !**************************************************************
    ! Create a separate ObjectItem for each handle

    subroutine createFpdbObjectHandle (rfpdbObjectItem, ihandle)

      ! The object item that represents the handle
      type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem

      ! The handle number
      integer, intent(in) :: ihandle


      ! local variables
      type(t_fpdbDataItem), pointer :: p_fpdbDataItem

      ! Check if storage block has UUID; otherwise create one
      if (uuid_isNil(rfpdbObjectItem%ruuid)) then
        call uuid_createUUID(4, rfpdbObjectItem%ruuid)
      end if
      
      ! Set the name and type of the object
      rfpdbObjectItem%sname = 'handle_'//trim(sys_siL(ihandle,11))
      rfpdbObjectItem%stype = 't_storageNode'

      ! Allocate the array of data items:
      allocate(rfpdbObjectItem%p_RfpdbDataItem(6))

      ! Fill the array of data items: ihandle
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
      p_fpdbDataItem%ctype    = FPDB_INT
      p_fpdbDataItem%sname    = 'ihandle'
      p_fpdbDataItem%iinteger = ihandle

      ! Fill the array of data items: idatatype
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
      p_fpdbDataItem%ctype    = FPDB_INT
      p_fpdbDataItem%sname    = 'idatatype'
      p_fpdbDataItem%iinteger = p_rheap%p_Rdescriptors(ihandle)%idatatype

      ! Fill the array of data items: idimension
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
      p_fpdbDataItem%ctype    = FPDB_INT
      p_fpdbDataItem%sname    = 'idimension'
      p_fpdbDataItem%iinteger = p_rheap%p_Rdescriptors(ihandle)%idimension

      ! Fill the array of data items: sname
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
      p_fpdbDataItem%ctype = FPDB_CHAR
      p_fpdbDataItem%sname = 'sname'
      p_fpdbDataItem%schar = p_rheap%p_Rdescriptors(ihandle)%sname

      ! Fill the array of data items: imemBytes
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
      p_fpdbDataItem%ctype  = FPDB_INT64
      p_fpdbDataItem%sname  = 'imemBytes'
      p_fpdbDataItem%iint64 = p_rheap%p_Rdescriptors(ihandle)%imemBytes

      
      ! Fill the array of data items: data
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
      p_fpdbDataItem%sname = 'data'

      select case(p_rheap%p_Rdescriptors(ihandle)%idimension)
      case (1)
        select case(p_rheap%p_Rdescriptors(ihandle)%idatatype)
        case (ST_SINGLE)
          p_fpdbDataItem%ctype       =  FPDB_SINGLE1D
          p_fpdbDataItem%p_Fsingle1D => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle1D

        case (ST_DOUBLE)
          p_fpdbDataItem%ctype       =  FPDB_DOUBLE1D
          p_fpdbDataItem%p_Ddouble1D => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble1D

        case (ST_QUAD)
          p_fpdbDataItem%ctype     =  FPDB_QUAD1D
          p_fpdbDataItem%p_Qquad1D => p_rheap%p_Rdescriptors(ihandle)%p_Qquad1D

        case (ST_INT)
          p_fpdbDataItem%ctype        =  FPDB_INT1D
          p_fpdbDataItem%p_Iinteger1D => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger1D

        case (ST_INT8)
          p_fpdbDataItem%ctype      =  FPDB_INT8_1D
          p_fpdbDataItem%p_Iint8_1D => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_1D

        case (ST_INT16)
          p_fpdbDataItem%ctype       =  FPDB_INT16_1D
          p_fpdbDataItem%p_Iint16_1D => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_1D

        case (ST_INT32)
          p_fpdbDataItem%ctype       =  FPDB_INT32_1D
          p_fpdbDataItem%p_Iint32_1D => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_1D

        case (ST_INT64)
          p_fpdbDataItem%ctype       =  FPDB_INT64_1D
          p_fpdbDataItem%p_Iint64_1D => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_1D

        case (ST_LOGICAL)
          p_fpdbDataItem%ctype        =  FPDB_LOGICAL1D
          p_fpdbDataItem%p_Blogical1D => p_rheap%p_Rdescriptors(ihandle)%p_Blogical1D

        case (ST_CHAR)
          p_fpdbDataItem%ctype        =  FPDB_CHAR1D
          p_fpdbDataItem%p_Schar1D => p_rheap%p_Rdescriptors(ihandle)%p_Schar1D

        case default
          call output_line ('Invalid data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'createFpdbObjectHandle')
          call sys_halt()
        end select

      case (2)
        select case(p_rheap%p_Rdescriptors(ihandle)%idatatype)
        case (ST_SINGLE)
          p_fpdbDataItem%ctype       =  FPDB_SINGLE2D
          p_fpdbDataItem%p_Fsingle2D => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle2D

        case (ST_DOUBLE)
          p_fpdbDataItem%ctype       =  FPDB_DOUBLE2D
          p_fpdbDataItem%p_Ddouble2D => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble2D

        case (ST_QUAD)
          p_fpdbDataItem%ctype     =  FPDB_QUAD2D
          p_fpdbDataItem%p_Qquad2D => p_rheap%p_Rdescriptors(ihandle)%p_Qquad2D

        case (ST_INT)
          p_fpdbDataItem%ctype        =  FPDB_INT2D
          p_fpdbDataItem%p_Iinteger2D => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger2D

        case (ST_INT8)
          p_fpdbDataItem%ctype      =  FPDB_INT8_2D
          p_fpdbDataItem%p_Iint8_2D => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_2D

        case (ST_INT16)
          p_fpdbDataItem%ctype       =  FPDB_INT16_2D
          p_fpdbDataItem%p_Iint16_2D => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_2D

        case (ST_INT32)
          p_fpdbDataItem%ctype       =  FPDB_INT32_2D
          p_fpdbDataItem%p_Iint32_2D => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_2D

        case (ST_INT64)
          p_fpdbDataItem%ctype       =  FPDB_INT64_2D
          p_fpdbDataItem%p_Iint64_2D => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_2D

        case (ST_LOGICAL)
          p_fpdbDataItem%ctype        =  FPDB_LOGICAL2D
          p_fpdbDataItem%p_Blogical2D => p_rheap%p_Rdescriptors(ihandle)%p_Blogical2D

        case (ST_CHAR)
          p_fpdbDataItem%ctype        =  FPDB_CHAR2D
          p_fpdbDataItem%p_Schar2D => p_rheap%p_Rdescriptors(ihandle)%p_Schar2D

        case default
          call output_line ('Invalid data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'createFpdbObjectHandle')
          call sys_halt()
        end select

      case (3)
          select case(p_rheap%p_Rdescriptors(ihandle)%idatatype)
        case (ST_SINGLE)
          p_fpdbDataItem%ctype       =  FPDB_SINGLE3D
          p_fpdbDataItem%p_Fsingle3D => p_rheap%p_Rdescriptors(ihandle)%p_Fsingle3D

        case (ST_DOUBLE)
          p_fpdbDataItem%ctype       =  FPDB_DOUBLE3D
          p_fpdbDataItem%p_Ddouble3D => p_rheap%p_Rdescriptors(ihandle)%p_Ddouble3D

        case (ST_QUAD)
          p_fpdbDataItem%ctype     =  FPDB_QUAD3D
          p_fpdbDataItem%p_Qquad3D => p_rheap%p_Rdescriptors(ihandle)%p_Qquad3D

        case (ST_INT)
          p_fpdbDataItem%ctype        =  FPDB_INT3D
          p_fpdbDataItem%p_Iinteger3D => p_rheap%p_Rdescriptors(ihandle)%p_Iinteger3D

        case (ST_INT8)
          p_fpdbDataItem%ctype      =  FPDB_INT8_3D
          p_fpdbDataItem%p_Iint8_3D => p_rheap%p_Rdescriptors(ihandle)%p_Iint8_3D

        case (ST_INT16)
          p_fpdbDataItem%ctype       =  FPDB_INT16_3D
          p_fpdbDataItem%p_Iint16_3D => p_rheap%p_Rdescriptors(ihandle)%p_Iint16_3D

        case (ST_INT32)
          p_fpdbDataItem%ctype       =  FPDB_INT32_3D
          p_fpdbDataItem%p_Iint32_3D => p_rheap%p_Rdescriptors(ihandle)%p_Iint32_3D

        case (ST_INT64)
          p_fpdbDataItem%ctype       =  FPDB_INT64_3D
          p_fpdbDataItem%p_Iint64_3D => p_rheap%p_Rdescriptors(ihandle)%p_Iint64_3D

        case (ST_LOGICAL)
          p_fpdbDataItem%ctype        =  FPDB_LOGICAL3D
          p_fpdbDataItem%p_Blogical3D => p_rheap%p_Rdescriptors(ihandle)%p_Blogical3D

        case (ST_CHAR)
          p_fpdbDataItem%ctype        =  FPDB_CHAR3D
          p_fpdbDataItem%p_Schar3D => p_rheap%p_Rdescriptors(ihandle)%p_Schar3D

        case default
          call output_line ('Invalid data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'createFpdbObjectHandle')
          call sys_halt()
        end select

      case default
        call output_line ('Invalid dimension!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'createFpdbObjectHandle')
        call sys_halt()
      end select
      
    end subroutine createFpdbObjectHandle
    
  end subroutine storage_createFpdbObject

!************************************************************************

!<subroutine>

  subroutine storage_restoreFpdbObject (rfpdbObjectItem, rheap)

!<description>
    ! This subroutine restores the heap structure from the abstract ObjectItem
!</description>

!<input>
    ! The object item that is created
    type(t_fpdbObjectItem), intent(in) :: rfpdbObjectItem
!</input>

!<inputoutput>
    ! OPTIONAL: local heap structure to initialise.
    ! If not given, the global heap is initialised.
    type(t_storageBlock), intent(inout), target, optional :: rheap
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_storageBlock), pointer :: p_rheap
    type(t_fpdbDataItem), pointer :: p_fpdbDataItem
    type(t_fpdbObjectItem), pointer :: p_fpdbObjectItem

    integer :: i

    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    ! Check if ObjectItem has correct type
    if (trim(rfpdbObjectItem%stype) .ne. 't_storageBlock') then
      call output_line ('Invalid object type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    end if

    ! Check if DataItems are associated
    if (.not.associated(rfpdbObjectItem%p_RfpdbDataItem)) then
      call output_line ('Missing data!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    end if

    ! Check if DataItems have correct size
    if (size(rfpdbObjectItem%p_RfpdbDataItem) .lt. 9) then
      call output_line ('Invalid data!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    end if

    ! Restore the data from the DataItem: p_inextFreeHandle
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'p_inextFreeHandle')) then
      call output_line ('Invalid data: p_inextFreeHandle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      p_rheap%p_inextFreeHandle = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: p_ilastFreeHandle
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'p_ilastFreeHandle')) then
      call output_line ('Invalid data: p_ilastFreeHandle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      p_rheap%p_ilastFreeHandle = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: ihandlesInUse
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'ihandlesInUse')) then
      call output_line ('Invalid data: ihandlesInUse!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      p_rheap%ihandlesInUse = p_fpdbDataItem%iinteger
    end if
    
    ! Restore the data from the DataItem: nhandlesTotal
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'nhandlesTotal')) then
      call output_line ('Invalid data: nhandlesTotal!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      p_rheap%nhandlesTotal = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: ihandlesDelta
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'ihandlesDelta')) then
      call output_line ('Invalid data: ihandlesDelta!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      p_rheap%ihandlesDelta = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: itotalMem
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT64) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'itotalMem')) then
      call output_line ('Invalid data: itotalMem!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      p_rheap%itotalMem = p_fpdbDataItem%iint64
    end if

    ! Restore the data from the DataItem: nhandlesInUseMax
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'nhandlesInUseMax')) then
      call output_line ('Invalid data: nhandlesInUseMax!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      p_rheap%nhandlesInUseMax = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: itotalMemMax
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT64) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'itotalMemMax')) then
      call output_line ('Invalid data: itotalMemMax!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      p_rheap%itotalMemMax = p_fpdbDataItem%iint64
    end if

    ! Restore the data from the DataItem: p_IfreeHandles
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT1D) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'p_IfreeHandles')) then
      call output_line ('Invalid data: p_IfreeHandles!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
      call sys_halt()
    else
      allocate(p_rheap%p_IfreeHandles(&
               p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
      call fpdb_getdata_int1d(p_fpdbDataItem, p_rheap%p_IfreeHandles)
    end if
    
    ! Initialise descriptors
    allocate(p_rheap%p_Rdescriptors(size(p_rheap%p_IfreeHandles)))

    ! Loop over remaining data and initialise memory nodes
    do i = 10, size(rfpdbObjectItem%p_RfpdbDataItem)
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(i)

      if (p_fpdbDataItem%ctype .ne. FPDB_OBJECT) then
        call output_line ('Datatype is no handle!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_restoreFpdbObject')
        call sys_halt()
      else
        p_fpdbObjectItem => p_fpdbDataItem%p_fpdbObjectItem
        call restoreFpdbObjectHandle(p_fpdbObjectItem)
      end if
    end do

  contains
    
    subroutine restoreFpdbObjectHandle(rfpdbObjectItem)

      ! The object item that represents the handle
      type(t_fpdbObjectItem), intent(in) :: rfpdbObjectItem

      
      ! local variables
      type(t_fpdbDataItem), pointer :: p_fpdbDataItem
      type(t_storageNode), pointer :: p_rnode
      integer :: ihandle


      ! Check if ObjectItem has correct type
      if (trim(rfpdbObjectItem%stype) .ne. 't_storageNode') then
        call output_line ('Invalid object type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      end if

      ! Check if DataItems are associated
      if (.not.associated(rfpdbObjectItem%p_RfpdbDataItem)) then
        call output_line ('Missing data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      end if

      ! Check if DataItems have correct size
      if (size(rfpdbObjectItem%p_RfpdbDataItem) .ne. 6) then
        call output_line ('Invalid data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      end if
      
      ! Restore the data from the DataItem: ihandle
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
      if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
          (trim(p_fpdbDataItem%sname) .ne. 'ihandle')) then
        call output_line ('Invalid data: ihandle!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      else
        ihandle =  p_fpdbDataItem%iinteger
        p_rnode => p_rheap%p_Rdescriptors(ihandle)
      end if

      ! Restore the data from the DataItem: idatatype
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
      if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
          (trim(p_fpdbDataItem%sname) .ne. 'idatatype')) then
        call output_line ('Invalid data: idatatype!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      else
        p_rnode%idatatype = p_fpdbDataItem%iinteger
      end if
      
      ! Restore the data from the DataItem: idimension
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
      if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
          (trim(p_fpdbDataItem%sname) .ne. 'idimension')) then
        call output_line ('Invalid data: idimension!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      else
        p_rnode%idimension = p_fpdbDataItem%iinteger
      end if

      ! Restore the data from the DataItem: sname
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
      if ((p_fpdbDataItem%ctype .ne. FPDB_CHAR) .or.&
          (trim(p_fpdbDataItem%sname) .ne. 'sname')) then
        call output_line ('Invalid data: sname!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      else
        p_rnode%sname = p_fpdbDataItem%schar
      end if

      ! Restore the data from the DataItem: imemBytes
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
      if ((p_fpdbDataItem%ctype .ne. FPDB_INT64) .or.&
          (trim(p_fpdbDataItem%sname) .ne. 'imemBytes')) then
        call output_line ('Invalid data: imemBytes!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      else
        p_rnode%imemBytes = p_fpdbDataItem%iint64
      end if


      ! Restore the data from the DataItem: data
      p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
      if (trim(p_fpdbDataItem%sname) .ne. 'data') then
        call output_line ('Invalid data: data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      end if

      select case(p_rnode%idimension)
      case (1)
        select case(p_rnode%idatatype)
        case (ST_SINGLE)
          allocate(p_rnode%p_Fsingle1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_single1d(p_fpdbDataItem, p_rnode%p_Fsingle1D)

        case (ST_DOUBLE)
          allocate(p_rnode%p_Ddouble1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_double1d(p_fpdbDataItem, p_rnode%p_Ddouble1D)

        case (ST_QUAD)
          allocate(p_rnode%p_Qquad1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_quad1d(p_fpdbDataItem, p_rnode%p_Qquad1D)

        case (ST_INT)
          allocate(p_rnode%p_Iinteger1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_int1d(p_fpdbDataItem, p_rnode%p_Iinteger1D)

        case (ST_INT8)
          allocate(p_rnode%p_Iint8_1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_int8_1d(p_fpdbDataItem, p_rnode%p_Iint8_1D)

        case (ST_INT16)
          allocate(p_rnode%p_Iint16_1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_int16_1d(p_fpdbDataItem, p_rnode%p_Iint16_1D)

        case (ST_INT32)
          allocate(p_rnode%p_Iint32_1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_int32_1d(p_fpdbDataItem, p_rnode%p_Iint32_1D)

        case (ST_INT64)
          allocate(p_rnode%p_Iint64_1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_int64_1d(p_fpdbDataItem, p_rnode%p_Iint64_1D)

        case (ST_LOGICAL)
          allocate(p_rnode%p_Blogical1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_logical1d(p_fpdbDataItem, p_rnode%p_Blogical1D)

        case (ST_CHAR)
          allocate(p_rnode%p_Schar1D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1)))
          call fpdb_getdata_char1d(p_fpdbDataItem, p_rnode%p_Schar1D)

        case DEFAULT
          call output_line ('Invalid data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
          call sys_halt()
        end select

      case (2)
        select case(p_rnode%idatatype)
        case (ST_SINGLE)
          allocate(p_rnode%p_Fsingle2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_single2d(p_fpdbDataItem, p_rnode%p_Fsingle2D)

        case (ST_DOUBLE)
          allocate(p_rnode%p_Ddouble2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_double2d(p_fpdbDataItem, p_rnode%p_Ddouble2D)

        case (ST_QUAD)
          allocate(p_rnode%p_Qquad2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_quad2d(p_fpdbDataItem, p_rnode%p_Qquad2D)

        case (ST_INT)
          allocate(p_rnode%p_Iinteger2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_int2d(p_fpdbDataItem, p_rnode%p_Iinteger2D)

        case (ST_INT8)
          allocate(p_rnode%p_Iint8_2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_int8_2d(p_fpdbDataItem, p_rnode%p_Iint8_2D)

        case (ST_INT16)
          allocate(p_rnode%p_Iint16_2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_int16_2d(p_fpdbDataItem, p_rnode%p_Iint16_2D)

        case (ST_INT32)
          allocate(p_rnode%p_Iint32_2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_int32_2d(p_fpdbDataItem, p_rnode%p_Iint32_2D)

        case (ST_INT64)
          allocate(p_rnode%p_Iint64_2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_int64_2d(p_fpdbDataItem, p_rnode%p_Iint64_2D)

        case (ST_LOGICAL)
          allocate(p_rnode%p_Blogical2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_logical2d(p_fpdbDataItem, p_rnode%p_Blogical2D)

        case (ST_CHAR)
          allocate(p_rnode%p_Schar2D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2)))
          call fpdb_getdata_char2d(p_fpdbDataItem, p_rnode%p_Schar2D)

        case default
          call output_line ('Invalid data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
          call sys_halt()
        end select

      case (3)
        select case(p_rnode%idatatype)
        case (ST_SINGLE)
          allocate(p_rnode%p_Fsingle3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_single3d(p_fpdbDataItem, p_rnode%p_Fsingle3D)

        case (ST_DOUBLE)
          allocate(p_rnode%p_Ddouble3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_double3d(p_fpdbDataItem, p_rnode%p_Ddouble3D)

        case (ST_QUAD)
          allocate(p_rnode%p_Qquad3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_quad3d(p_fpdbDataItem, p_rnode%p_Qquad3D)

        case (ST_INT)
          allocate(p_rnode%p_Iinteger3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_int3d(p_fpdbDataItem, p_rnode%p_Iinteger3D)

        case (ST_INT8)
          allocate(p_rnode%p_Iint8_3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_int8_3d(p_fpdbDataItem, p_rnode%p_Iint8_3D)

        case (ST_INT16)
          allocate(p_rnode%p_Iint16_3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_int16_3d(p_fpdbDataItem, p_rnode%p_Iint16_3D)

        case (ST_INT32)
          allocate(p_rnode%p_Iint32_3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_int32_3d(p_fpdbDataItem, p_rnode%p_Iint32_3D)

        case (ST_INT64)
          allocate(p_rnode%p_Iint64_3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_int64_3d(p_fpdbDataItem, p_rnode%p_Iint64_3D)

        case (ST_LOGICAL)
          allocate(p_rnode%p_Blogical3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_logical3d(p_fpdbDataItem, p_rnode%p_Blogical3D)

        case (ST_CHAR)
          allocate(p_rnode%p_Schar3D(&
                   p_fpdbDataItem%Ilbounds(1):p_fpdbDataItem%Iubounds(1),&
                   p_fpdbDataItem%Ilbounds(2):p_fpdbDataItem%Iubounds(2),&
                   p_fpdbDataItem%Ilbounds(3):p_fpdbDataItem%Iubounds(3)))
          call fpdb_getdata_char3d(p_fpdbDataItem, p_rnode%p_Schar3D)

        case default
          call output_line ('Invalid data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
          call sys_halt()
        end select

      case default
        call output_line ('Invalid dimension!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectHandle')
        call sys_halt()
      end select
      
    end subroutine restoreFpdbObjectHandle

  end subroutine storage_restoreFpdbObject
  
!************************************************************************

!<subroutine>

  subroutine storage_setdatatype (ihandle, ctype, rheap)

!<description>
  ! Converts an array to another data type.
  ! Note that this usually involves a reallocation of the
  ! memory and probably a loss in accuracy.
  ! Not all conversions are possible.
  ! Converting to/from a string is not possible.
!</description>

!<input>
  ! Handle of the memory block to be releases
  integer, intent(in) :: ihandle

  ! Target Datatype of the array identified by ihandle.
  integer, intent(in) :: ctype

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</input>

!</subroutine>

  ! local variables.
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode
  integer :: inewhandle,isize1d
  integer, dimension(:), allocatable :: Isize
  type(t_storageNode) :: rtempnode

    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if

    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Handle invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_setdatatype')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    if ((ctype .eq. ST_CHAR) .or. (p_rnode%idataType .eq. ST_CHAR)) then
      call output_line ('Converting to/from a string is not allowed!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_setdatatype')
      call sys_halt()
    end if
    
    ! Probably nothing to do.
    if (ctype .eq. p_rnode%idataType) return

    ! Copy our storage node to an empty storage node.
    select case (p_rnode%idimension)
    case (1)
      call storage_getsize (ihandle,isize1d,rheap)
      call storage_new ('storage_setdatatype', p_rnode%sname, isize1d, &
          ctype, inewhandle, ST_NEWBLOCK_NOINIT, rheap)
    case default
      allocate(Isize(p_rnode%idimension))
      call storage_getsize (ihandle,Isize,rheap)
      call storage_new ('storage_setdatatype', p_rnode%sname, Isize, &
          ctype, inewhandle, ST_NEWBLOCK_NOINIT, rheap)
      deallocate(Isize)
    end select
    
    ! Copy the data. This will do the data conversion.
    call storage_copy (ihandle,inewhandle)
    
    ! Exchange the two storage nodes.
    rtempnode = p_rheap%p_Rdescriptors(ihandle)
    p_rheap%p_Rdescriptors(ihandle) = p_rheap%p_Rdescriptors(inewhandle)
    p_rheap%p_Rdescriptors(inewhandle) = rtempnode
    
    ! Release th old node which has now the new handle as number.
    call storage_free (inewhandle)
  
  end subroutine storage_setdatatype

!************************************************************************

!<function>

  function storage_getblocktype (stype) result (ctype)

!<description>
  ! This function tries to determine the corresponding data type
  !  ST_XXX from the string stype and return ST_NOHANDLE otherwise.
!</description>

!<input>
  character(LEN=*), intent(in) :: stype
!</input>

!<result>
  integer :: ctype
!</result>

!</function>

    ! local variable
    character(len=(len(stype))) :: sdatatype

    call sys_tolower(stype, sdatatype)

    ! Determine data type
    select case((trim(adjustl(sdatatype))))
    case ('single','real(sp)')
      ctype = ST_SINGLE
    case ('double','real(dp)')
      ctype = ST_DOUBLE
    case ('quad','real(qp)')
      ctype = ST_QUAD
    case ('integer')
      ctype = ST_INT
    case ('int8','integer(i8)')
      ctype = ST_INT8
    case ('int16','integer(i16)')
      ctype = ST_INT16
    case ('int32','integer(i32)')
      ctype = ST_INT32
    case ('int64','integer(i64)')
      ctype = ST_INT64
    case ('logical')
      ctype = ST_LOGICAL
    case ('char','character')
      ctype = ST_CHAR
    case default
      ctype = ST_NOHANDLE
    end select
  end function storage_getblocktype

!************************************************************************

!<subroutine>

  subroutine storage_allocMemoryOnDevice (ihandle, rheap)

!<description>
  ! This routine allocates a memory block in the device memory
  ! associated with the given handle of the heap. The handle
  ! must already be associated with some memory block.
!</description>

!<input>
  ! Handle of the memory block
  integer, intent(in) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure.
  ! If not given, the global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</inputoutput>

!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

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
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_allocMemoryOnDevice')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Check if memory pointer is already in use
    if (storage_isAssociated(p_rnode%cdeviceMemPtr))&
        call coproc_freeMemoryOnDevice(p_rnode%cdeviceMemPtr)

    ! Allocate memory on device with correct size
    call coproc_newMemoryOnDevice(p_rnode%cdeviceMemPtr, p_rnode%imemBytes)

#else

    call output_line ('Application must be compiled with coprocessor support enabled!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_allocMemoryOnDevice')
    call sys_halt()

#endif

  end subroutine storage_allocMemoryOnDevice

!************************************************************************

!<subroutine>

  subroutine storage_deallocMemoryOnDevice (ihandle, rheap)

!<description>
  ! This routine deallocates a memory block in the device memory
  ! associated with the given handle of the heap. The handle must
  ! already be associated with some memory block.
!</description>

!<input>
  ! Handle of the memory block
  integer, intent(in) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure.
  ! If not given, the global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</inputoutput>

!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

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
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_deallocMemoryOnDevice')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Check if memory pointer is already in use
    if (storage_isAssociated(p_rnode%cdeviceMemPtr))&
        call coproc_freeMemoryOnDevice(p_rnode%cdeviceMemPtr, p_rnode%imemBytes)

#else

    call output_line ('Application must be compiled with coprocessor support enabled!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_deallocMemoryOnDevice')
    call sys_halt()

#endif

  end subroutine storage_deallocMemoryOnDevice

!************************************************************************

!<subroutine>


  recursive subroutine storage_syncMemoryHostDevice (ihandle, csyncBlock,&
      btransposeMemoryOrder, istream, rheap)

!<description>
  ! This routine synchronises the handle of the heap between host
  ! memory and device memory.
!</description>

!<input>
  ! Handle of the memory block
  integer, intent(in) :: ihandle

  ! Synchronisation identifier (ST_SYNCBLOCK_COPY_H2D,
  ! ST_SYNCBLOCK_COPY_D2H, ST_SYNCBLOCK_ACCUMULATE_D2H,
  ! ST_SYNCBLOCK_ACCUMULATE_H2D). Specifies how to
  ! synchronise the data in host and device memory.
  integer, intent(in) :: csyncBlock

  ! OPTIONAL: if .TRUE. then the memory order is transposed before
  ! transfering it from host to device memory and vice versa
  logical, intent(in), optional :: btransposeMemoryOrder

  ! OPTIONAL: stream for asynchronious transfer.
  ! If istream is present and if asynchroneous transfer is supported
  ! then all memory transfers are carried out asynchroneously
  integer(I64), intent(in), optional :: istream
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</inputoutput>

!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

  ! Pointer to the heap
  type(t_storageBlock), pointer :: p_rheap
  type(t_storageNode), pointer :: p_rnode, p_rnodeTmp
  type(C_PTR) :: cdeviceMemPtr
  integer :: ihandleTmp
  logical :: btranspose
  integer(I64) :: istreamTmp
  
    ! Do we have to transpose the memory order
    btranspose = .false.
    if (present(btransposeMemoryOrder)) btranspose = btransposeMemoryOrder
  
    ! Get the heap to use - local or global one.
    if(present(rheap)) then
      p_rheap => rheap
    else
      p_rheap => rbase
    end if
    
    if (ihandle .le. ST_NOHANDLE) then
      call output_line ('Wrong handle!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    select case(csyncBlock)
    case (ST_SYNCBLOCK_COPY_H2D)

      ! Copy memory block associated with handle ihandle from host
      ! memory to device memory; if memory is not allocated on device
      ! then a new memory block is first allocated on the device

      ! Check if memory on device is allocated
      if (.not.storage_isAssociated(p_rnode%cdeviceMemPtr))&
          call coproc_newMemoryOnDevice(p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
      
#ifdef USE_C_PTR_STORAGE
      ! Check if asynchroneous transfer is applicable
      if (present(istream) .and. .not.(btranspose)) then
        call coproc_memcpyHostToDeviceAsync(p_rnode%chostMemPtr,&
            p_rnode%cdeviceMemPtr, p_rnode%imemBytes, istream)
        
        ! That`s it
        return
      end if
#endif
      
      ! What dimension are we?
      select case(p_rnode%idimension)
        
      case (1)
        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          call  coproc_memcpyHostToDevice(p_rnode%p_Fsingle1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_DOUBLE)
          call coproc_memcpyHostToDevice(p_rnode%p_Ddouble1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_QUAD)
          call coproc_memcpyHostToDevice(p_rnode%p_Qquad1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_INT)
          call coproc_memcpyHostToDevice(p_rnode%p_Iinteger1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_INT8)
          call coproc_memcpyHostToDevice(p_rnode%p_Iint8_1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_INT16)
          call coproc_memcpyHostToDevice(p_rnode%p_Iint16_1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_INT32)
          call coproc_memcpyHostToDevice(p_rnode%p_Iint32_1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_INT64)
          call coproc_memcpyHostToDevice(p_rnode%p_Iint64_1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_LOGICAL)
          call coproc_memcpyHostToDevice(p_rnode%p_Blogical1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case (ST_CHAR)
          call coproc_memcpyHostToDevice(p_rnode%p_Schar1D,&
              p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
        case default
          call output_line ('Unsupported data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select

      case (2)
        ! Do we have to transpose the memory order?
        if (btranspose) then
          ! TODO: Once an in-place transpose has been implemented on
          ! the coprocessor device, we do not need this temporal
          ! storage but can copy the original data to the device and
          ! perform the transpose operation directly on the device. So
          ! far, we have to work with a temporal array since the data
          ! may be needed untranspoed in host memory.

          ! Create a temporal working array ...
          ihandleTmp = ST_NOHANDLE
          call storage_new('storage_syncMemoryHostDevice', 'ihandleTmp',&
                           ihandle, ihandleTmp, p_rheap)
          ! ... and get the associated storage node
          p_rnodeTmp => p_rheap%p_Rdescriptors(ihandleTmp)
        end if

        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          if (btranspose) then
            call transposeOutOfPlace_Single2D(p_rnode%p_Fsingle2D, p_rnodeTmp%p_Fsingle2D,&
                size(p_rnodeTmp%p_Fsingle2D,1), size(p_rnodeTmp%p_Fsingle2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Fsingle2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Fsingle2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_DOUBLE)
          if (btranspose) then
            call transposeOutOfPlace_Double2D(p_rnode%p_Ddouble2D, p_rnodeTmp%p_Ddouble2D,&
                size(p_rnodeTmp%p_Ddouble2D,1), size(p_rnodeTmp%p_Ddouble2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Ddouble2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Ddouble2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_QUAD)
          if (btranspose) then
            call transposeOutOfPlace_Quad2D(p_rnode%p_Qquad2D, p_rnodeTmp%p_Qquad2D,&
                size(p_rnodeTmp%p_Qquad2D,1), size(p_rnodeTmp%p_Qquad2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Qquad2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Qquad2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT)
          if (btranspose) then
            call transposeOutOfPlace_Integer2D(p_rnode%p_Iinteger2D, p_rnodeTmp%p_Iinteger2D,&
                size(p_rnodeTmp%p_Iinteger2D,1), size(p_rnodeTmp%p_Iinteger2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iinteger2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iinteger2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT8)
          if (btranspose) then
            call transposeOutOfPlace_Int8_2D(p_rnode%p_Iint8_2D, p_rnodeTmp%p_Iint8_2D,&
                size(p_rnodeTmp%p_Iint8_2D,1), size(p_rnodeTmp%p_Iint8_2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iint8_2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iint8_2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT16)
          if (btranspose) then
            call transposeOutOfPlace_Int16_2D(p_rnode%p_Iint16_2D, p_rnodeTmp%p_Iint16_2D,&
                size(p_rnodeTmp%p_Iint16_2D,1), size(p_rnodeTmp%p_Iint16_2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iint16_2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iint16_2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT32)
          if (btranspose) then
            call transposeOutOfPlace_Int32_2D(p_rnode%p_Iint32_2D, p_rnodeTmp%p_Iint32_2D,&
                size(p_rnodeTmp%p_Iint32_2D,1), size(p_rnodeTmp%p_Iint32_2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iint32_2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iint32_2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT64)
          if (btranspose) then
            call transposeOutOfPlace_Int64_2D(p_rnode%p_Iint64_2D, p_rnodeTmp%p_Iint64_2D,&
                size(p_rnodeTmp%p_Iint64_2D,1), size(p_rnodeTmp%p_Iint64_2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iint64_2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iint64_2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_LOGICAL)
          if (btranspose) then
            call transposeOutOfPlace_Logical2D(p_rnode%p_Blogical2D, p_rnodeTmp%p_Blogical2D,&
                size(p_rnodeTmp%p_Blogical2D,1), size(p_rnodeTmp%p_Blogical2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Blogical2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Blogical2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_CHAR)
          if (btranspose) then
            call transposeOutOfPlace_Char2D(p_rnode%p_SChar2D, p_rnodeTmp%p_SChar2D,&
                size(p_rnodeTmp%p_SChar2D,1), size(p_rnodeTmp%p_SChar2D,2))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_SChar2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Schar2D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case default
          call output_line ('Unsupported data typee!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select

        ! Free temporal working array (if created)
        if (btranspose) call storage_free(ihandleTmp)

      case (3)
        ! Do we have to transpose the memory order?
        if (btranspose) then
          ! TODO: Once an in-place transpose has been implemented on
          ! the coprocessor device, we do not need this temporal
          ! storage but can copy the original data to the device and
          ! perform the transpose operation directly on the device. So
          ! far, we habe to work with a temporal array since the data
          ! may be needed untranspoed in host memory.

          ! Create a temporal working array ...
          ihandleTmp = ST_NOHANDLE
          call storage_new('storage_syncMemoryHostDevice', 'ihandleTmp',&
                           ihandle, ihandleTmp, p_rheap)
          ! ... and get the associated storage node
          p_rnodeTmp => p_rheap%p_Rdescriptors(ihandleTmp)
        end if

        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          if (btranspose) then
            call transposeOutOfPlace_Single3D(p_rnode%p_Fsingle3D, p_rnodeTmp%p_Fsingle3D,&
                size(p_rnodeTmp%p_Fsingle3D,1), size(p_rnodeTmp%p_Fsingle3D,2),&
                size(p_rnodeTmp%p_Fsingle3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Fsingle3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Fsingle3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_DOUBLE)
          if (btranspose) then
            call transposeOutOfPlace_Double3D(p_rnode%p_Ddouble3D, p_rnodeTmp%p_Ddouble3D,&
                size(p_rnodeTmp%p_Ddouble3D,1), size(p_rnodeTmp%p_Ddouble3D,2),&
                size(p_rnodeTmp%p_Ddouble3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Ddouble3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Ddouble3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_QUAD)
          if (btranspose) then
            call transposeOutOfPlace_Quad3D(p_rnode%p_Qquad3D, p_rnodeTmp%p_Qquad3D,&
                size(p_rnodeTmp%p_Qquad3D,1), size(p_rnodeTmp%p_Qquad3D,2),&
                size(p_rnodeTmp%p_Qquad3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Qquad3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Qquad3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT)
          if (btranspose) then
            call transposeOutOfPlace_Integer3D(p_rnode%p_Iinteger3D, p_rnodeTmp%p_Iinteger3D,&
                size(p_rnodeTmp%p_Iinteger3D,1), size(p_rnodeTmp%p_Iinteger3D,2),&
                size(p_rnodeTmp%p_Iinteger3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iinteger3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iinteger3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT8)
          if (btranspose) then
            call transposeOutOfPlace_Int8_3D(p_rnode%p_Iint8_3D, p_rnodeTmp%p_Iint8_3D,&
                size(p_rnodeTmp%p_Iint8_3D,1), size(p_rnodeTmp%p_Iint8_3D,2),&
                size(p_rnodeTmp%p_Iint8_3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iint8_3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iint8_3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT16)
          if (btranspose) then
            call transposeOutOfPlace_Int16_3D(p_rnode%p_Iint16_3D, p_rnodeTmp%p_Iint16_3D,&
                size(p_rnodeTmp%p_Iint16_3D,1), size(p_rnodeTmp%p_Iint16_3D,2),&
                size(p_rnodeTmp%p_Iint16_3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iint16_3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iint16_3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT32)
          if (btranspose) then
            call transposeOutOfPlace_Int32_3D(p_rnode%p_Iint32_3D, p_rnodeTmp%p_Iint32_3D,&
                size(p_rnodeTmp%p_Iint32_3D,1), size(p_rnodeTmp%p_Iint32_3D,2),&
                size(p_rnodeTmp%p_Iint32_3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iint32_3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iint32_3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_INT64)
          if (btranspose) then
            call transposeOutOfPlace_Int64_3D(p_rnode%p_Iint64_3D, p_rnodeTmp%p_Iint64_3D,&
                size(p_rnodeTmp%p_Iint64_3D,1), size(p_rnodeTmp%p_Iint64_3D,2),&
                size(p_rnodeTmp%p_Iint64_3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Iint64_3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Iint64_3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_LOGICAL)
          if (btranspose) then
            call transposeOutOfPlace_Logical3D(p_rnode%p_Blogical3D, p_rnodeTmp%p_Blogical3D,&
                size(p_rnodeTmp%p_Blogical3D,1), size(p_rnodeTmp%p_Blogical3D,2),&
                size(p_rnodeTmp%p_Blogical3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_Blogical3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Blogical3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case (ST_CHAR)
          if (btranspose) then
            call transposeOutOfPlace_Char3D(p_rnode%p_SChar3D, p_rnodeTmp%p_SChar3D,&
                size(p_rnodeTmp%p_SChar3D,1), size(p_rnodeTmp%p_SChar3D,2),&
                size(p_rnodeTmp%p_SChar3D,3))
            call coproc_memcpyHostToDevice(p_rnodeTmp%p_SChar3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          else
            call coproc_memcpyHostToDevice(p_rnode%p_Schar3D,&
                p_rnode%cdeviceMemPtr, p_rnode%imemBytes)
          end if

        case default
          call output_line ('Unsupported data typee!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select

        ! Free temporal working array (if created)
        if (btranspose) call storage_free(ihandleTmp)

      case default
        call output_line ('Unsupported dimension!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
        call sys_halt()
      end select

    case (ST_SYNCBLOCK_COPY_D2H)

      ! Copy memory block associated with handle ihandle from device
      ! memory to host memory; if memory is not allocated on device
      ! then an error is thrown and the program terminates.

      ! Check if memory on device is allocated
      if (.not.storage_isAssociated(p_rnode%cdeviceMemPtr)) then
        call output_line ('Invalid memory address!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
        call sys_halt()
      end if

#ifdef USE_C_PTR_STORAGE
        ! Check if asynchroneous transfer is applicable
        if (present(istream) .and. .not.(btranspose)) then
          call coproc_memcpyDeviceToHostAsync(p_rnode%cdeviceMemPtr,&
              p_rnode%chostMemPtr, p_rnode%imemBytes, istream)
          
          ! That`s it
          return
        end if
#endif

      ! What dimension are we?
      select case(p_rnode%idimension)
        
      case (1)
        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Fsingle1D, p_rnode%imemBytes)
        case (ST_DOUBLE)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Ddouble1D, p_rnode%imemBytes)
        case (ST_QUAD)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Qquad1D, p_rnode%imemBytes)
        case (ST_INT)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iinteger1D, p_rnode%imemBytes)
        case (ST_INT8)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint8_1D, p_rnode%imemBytes)
        case (ST_INT16)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint16_1D, p_rnode%imemBytes)
        case (ST_INT32)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint32_1D, p_rnode%imemBytes)
        case (ST_INT64)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint64_1D, p_rnode%imemBytes)
        case (ST_LOGICAL)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Blogical1D, p_rnode%imemBytes)
        case (ST_CHAR)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Schar1D, p_rnode%imemBytes)
        case default
          call output_line ('Unsupported data typee!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select

      case (2)
        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Fsingle2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Single2D(p_rnode%p_Fsingle2D,&
              size(p_rnode%p_Fsingle2D,2), size(p_rnode%p_Fsingle2D,1))

        case (ST_DOUBLE)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Ddouble2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Double2D(p_rnode%p_Ddouble2D,&
              size(p_rnode%p_Ddouble2D,2), size(p_rnode%p_Ddouble2D,1))

        case (ST_QUAD)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Qquad2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Quad2D(p_rnode%p_Qquad2D,&
              size(p_rnode%p_Qquad2D,2), size(p_rnode%p_Qquad2D,1))

        case (ST_INT)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iinteger2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Integer2D(p_rnode%p_Iinteger2D,&
              size(p_rnode%p_Iinteger2D,2), size(p_rnode%p_Iinteger2D,1))

        case (ST_INT8)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint8_2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Int8_2D(p_rnode%p_Iint8_2D,&
              size(p_rnode%p_Iint8_2D,2), size(p_rnode%p_Iint8_2D,1))

        case (ST_INT16)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint16_2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Int16_2D(p_rnode%p_Iint16_2D,&
              size(p_rnode%p_Iint16_2D,2), size(p_rnode%p_Iint16_2D,1))

        case (ST_INT32)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint32_2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Int32_2D(p_rnode%p_Iint32_2D,&
              size(p_rnode%p_Iint32_2D,2), size(p_rnode%p_Iint32_2D,1))

        case (ST_INT64)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint64_2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Int64_2D(p_rnode%p_Iint64_2D,&
              size(p_rnode%p_Iint64_2D,2), size(p_rnode%p_Iint64_2D,1))

        case (ST_LOGICAL)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Blogical2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Logical2D(p_rnode%p_Blogical2D,&
              size(p_rnode%p_Blogical2D,2), size(p_rnode%p_Blogical2D,1))

        case (ST_CHAR)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Schar2D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Char2D(p_rnode%p_Schar2D,&
              size(p_rnode%p_Schar2D,2), size(p_rnode%p_Schar2D,1))

        case default
          call output_line ('Unsupported data typee!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select

      case (3)
        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Fsingle3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Single3D(p_rnode%p_Fsingle3D,&
              size(p_rnode%p_Fsingle3D,3), size(p_rnode%p_Fsingle3D,2),&
              size(p_rnode%p_Fsingle3D,1))

        case (ST_DOUBLE)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Ddouble3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Double3D(p_rnode%p_Ddouble3D,&
              size(p_rnode%p_Ddouble3D,3), size(p_rnode%p_Ddouble3D,2),&
              size(p_rnode%p_Ddouble3D,1))

        case (ST_QUAD)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Qquad3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Quad3D(p_rnode%p_Qquad3D,&
              size(p_rnode%p_Qquad3D,3), size(p_rnode%p_Qquad3D,2),&
              size(p_rnode%p_Qquad3D,1))

        case (ST_INT)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iinteger3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Integer3D(p_rnode%p_Iinteger3D,&
              size(p_rnode%p_Iinteger3D,3), size(p_rnode%p_Iinteger3D,2),&
              size(p_rnode%p_Iinteger3D,1))

        case (ST_INT8)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint8_3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Int8_3D(p_rnode%p_Iint8_3D,&
              size(p_rnode%p_Iint8_3D,3), size(p_rnode%p_Iint8_3D,2),&
              size(p_rnode%p_Iint8_3D,1))

        case (ST_INT16)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint16_3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Int16_3D(p_rnode%p_Iint16_3D,&
              size(p_rnode%p_Iint16_3D,3), size(p_rnode%p_Iint16_3D,2),&
              size(p_rnode%p_Iint16_3D,1))

        case (ST_INT32)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint32_3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Int32_3D(p_rnode%p_Iint32_3D,&
              size(p_rnode%p_Iint32_3D,3), size(p_rnode%p_Iint32_3D,2),&
              size(p_rnode%p_Iint32_3D,1))

        case (ST_INT64)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Iint64_3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Int64_3D(p_rnode%p_Iint64_3D,&
              size(p_rnode%p_Iint64_3D,3), size(p_rnode%p_Iint64_3D,2),&
              size(p_rnode%p_Iint64_3D,1))

        case (ST_LOGICAL)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Blogical3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Logical3D(p_rnode%p_Blogical3D,&
              size(p_rnode%p_Blogical3D,3), size(p_rnode%p_Blogical3D,2),&
              size(p_rnode%p_Blogical3D,1))

        case (ST_CHAR)
          call coproc_memcpyDeviceToHost(p_rnode%cdeviceMemPtr,&
              p_rnode%p_Schar3D, p_rnode%imemBytes)
          if (btranspose) call transposeInPlace_Char3D(p_rnode%p_Schar3D,&
              size(p_rnode%p_Schar3D,3), size(p_rnode%p_Schar3D,2),&
              size(p_rnode%p_Schar3D,1))

        case default
          call output_line ('Unsupported data typee!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select

      case default
        call output_line ('Unsupported dimension!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
        call sys_halt()
      end select

    case (ST_SYNCBLOCK_ACCUMULATE_H2D)

      ! Copy-add memory block associated with handle ihandle from host
      ! memory to device memory; if memory is not allocated on device
      ! then a new memory block is first allocated on the device

      ! Check if memory on device is allocated
      if (.not.storage_isAssociated(p_rnode%cdeviceMemPtr)) then
        ! Device memory has not been allocated so this is a simple copy
        call storage_syncMemoryHostDevice(ihandle, ST_SYNCBLOCK_COPY_H2D,&
            btranspose, istream, p_rheap)
      else
        ! This is tricky. We need to transfer the data from host memory
        ! into device memory and accumulate it into device memory. This
        ! cannot be done directly but some temporal storage is required.
      
        ! Make a backup of the memory pointer
        cdeviceMemPtr = p_rnode%cdeviceMemPtr
        
        ! Allocate memory on the device by hand
        call coproc_newMemoryOnDevice(p_rnode%cdeviceMemPtr, p_rnode%imemBytes)

        ! Copy content of host memory associated with handle ihandle
        ! into temporal memory block on device. Note that the memory
        ! address has been modified by hand, so that data is preserved.
        call storage_syncMemoryHostDevice(ihandle, ST_SYNCBLOCK_COPY_H2D,&
            btranspose, istream, p_rheap)

        ! We are back from recursion and both memory blocks are
        ! available in device memory. So we can add both memory blocks
        ! and store the result at the original memory address.

        istreamTmp = 0_I64
        if (present(istream)) istreamTmp = istream

        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          call coproc_combineSingleOnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_SINGLE_BYTES, istreamTmp)
        case (ST_DOUBLE)
          call coproc_combineDoubleOnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_DOUBLE_BYTES, istreamTmp)
        case (ST_QUAD)
          call coproc_combineQuadOnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_QUAD_BYTES, istreamTmp)
        case (ST_INT)
          call coproc_combineIntegerOnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_INT_BYTES, istreamTmp)
        case (ST_INT8)
          call coproc_combineInt8OnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_INT8_BYTES, istreamTmp)
        case (ST_INT16)
          call coproc_combineInt16OnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_INT16_BYTES, istreamTmp)
        case (ST_INT32)
          call coproc_combineInt32OnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_INT32_BYTES, istreamTmp)
        case (ST_INT64)
          call coproc_combineInt64OnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_INT64_BYTES, istreamTmp)
        case (ST_LOGICAL)
          call coproc_combineLogicalOnDevice(p_rnode%cdeviceMemPtr, cdeviceMemPtr,&
              cdeviceMemPtr, p_rnode%imemBytes/ST_LOGICAL_BYTES, istreamTmp)
        case default
          call output_line ('Unsupported data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
        end select

        ! Synchronise stream
        call coproc_synchroniseStream(istreamTmp)
      
        ! Release temporal memory block ...
        call coproc_freeMemoryOnDevice(p_rnode%cdeviceMemPtr)

        ! ... and restore backup of original memory address
        p_rnode%cdeviceMemPtr = cdeviceMemPtr
      end if

    case (ST_SYNCBLOCK_ACCUMULATE_D2H)
      
      ! Copy-add memory block associated with handle ihandle from device
      ! memory to host memory; if memory is not allocated on device
      ! then a new memory block is first allocated on the device

      ! Check if memory on device is allocated
      if (.not.storage_isAssociated(p_rnode%cdeviceMemPtr)) then
        call output_line ('Invalid memory address!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
        call sys_halt()
      end if

      ! This is tricky. We need to transfer the data from device
      ! memory into host memory and accumulate it into host memory.
      
      ! For this we need some temporal storage in host memory
      ihandleTmp = ST_NOHANDLE
      call storage_new('storage_syncMemoryHostDevice', 'ihandleTmp',&
                       ihandle, ihandleTmp, p_rheap)

      ! Transfer the memory address associated with handle ihandle
      ! to the temporal memory block associated with ihandleTmp
      p_rnodeTmp => p_rheap%p_Rdescriptors(ihandleTmp)
      p_rnodeTmp%cdeviceMemPtr = p_rnode%cdeviceMemPtr
      
      ! Copy content from device memory to host memory block
      ! associated with temporal handle ihandleTmp. Note that it does
      ! not make sense to allow for asynchroneous transfers here,
      ! since the memory transfer must be complete before the two
      ! memory blocks in host memory can be combined.
      call storage_syncMemoryHostDevice(ihandleTmp, ST_SYNCBLOCK_COPY_D2H,&
          btranspose, rheap=p_rheap)
      
      ! We are back from recursion and both memory blocks are
      ! available in host memory. So we can add both memory blocks
      ! and store the result at the original memory address.
            
      ! What dimension are we?
      select case(p_rnode%idimension)
        
      case (1)
        
        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          call lalg_vectorLinearCombSngl (p_rnodeTmp%p_Fsingle1D,&
              p_rnode%p_Fsingle1D, 1.0_SP, 1.0_SP)
        case (ST_DOUBLE)
          call lalg_vectorLinearCombDble (p_rnodeTmp%p_Ddouble1D,&
              p_rnode%p_Ddouble1D, 1.0_DP, 1.0_DP)
        case (ST_QUAD)
          call lalg_vectorLinearCombQuad (p_rnodeTmp%p_Qquad1D,&
              p_rnode%p_Qquad1D, 1.0_QP, 1.0_QP)
        case (ST_INT)
          p_rnode%p_Iinteger1D = p_rnode%p_Iinteger1D + p_rnodeTmp%p_Iinteger1D
        case (ST_INT8)
          p_rnode%p_Iint8_1D = p_rnode%p_Iint8_1D + p_rnodeTmp%p_Iint8_1D
        case (ST_INT16)
          p_rnode%p_Iint16_1D = p_rnode%p_Iint16_1D + p_rnodeTmp%p_Iint16_1D
        case (ST_INT32)
          p_rnode%p_Iint32_1D = p_rnode%p_Iint32_1D + p_rnodeTmp%p_Iint32_1D
        case (ST_INT64)
          p_rnode%p_Iint64_1D = p_rnode%p_Iint64_1D + p_rnodeTmp%p_Iint64_1D
        case (ST_LOGICAL)
          p_rnode%p_Blogical1D = p_rnode%p_Blogical1D .or. p_rnodeTmp%p_Blogical1D
        case default
          call output_line ('Unsupported data typee!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select
        
      case (2)
        
        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          call lalg_vectorLinearCombSngl2D (p_rnodeTmp%p_Fsingle2D,&
              p_rnode%p_Fsingle2D, 1.0_SP, 1.0_SP)
        case (ST_DOUBLE)
          call lalg_vectorLinearCombDble2D (p_rnodeTmp%p_Ddouble2D,&
              p_rnode%p_Ddouble2D, 1.0_DP, 1.0_DP)
        case (ST_QUAD)
          call lalg_vectorLinearCombQuad2D (p_rnodeTmp%p_Qquad2D,&
              p_rnode%p_Qquad2D, 1.0_QP, 1.0_QP)
        case (ST_INT)
          p_rnode%p_Iinteger2D = p_rnode%p_Iinteger2D + p_rnodeTmp%p_Iinteger2D
        case (ST_INT8)
          p_rnode%p_Iint8_2D = p_rnode%p_Iint8_2D + p_rnodeTmp%p_Iint8_2D
        case (ST_INT16)
          p_rnode%p_Iint16_2D = p_rnode%p_Iint16_2D + p_rnodeTmp%p_Iint16_2D
        case (ST_INT32)
          p_rnode%p_Iint32_2D = p_rnode%p_Iint32_2D + p_rnodeTmp%p_Iint32_2D
        case (ST_INT64)
          p_rnode%p_Iint64_2D = p_rnode%p_Iint64_2D + p_rnodeTmp%p_Iint64_2D
        case (ST_LOGICAL)
          p_rnode%p_Blogical2D = p_rnode%p_Blogical2D .or. p_rnodeTmp%p_Blogical2D
        case default
          call output_line ('Unsupported data typee!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select
        
      case (3)
        
        ! What data type are we?
        select case(p_rnode%idataType)
        case (ST_SINGLE)
          call lalg_vectorLinearCombSngl3D (p_rnodeTmp%p_Fsingle3D,&
              p_rnode%p_Fsingle3D, 1.0_SP, 1.0_SP)
        case (ST_DOUBLE)
          call lalg_vectorLinearCombDble3D (p_rnodeTmp%p_Ddouble3D,&
              p_rnode%p_Ddouble3D, 1.0_DP, 1.0_DP)
        case (ST_QUAD)
          call lalg_vectorLinearCombQuad3D (p_rnodeTmp%p_Qquad3D,&
              p_rnode%p_Qquad3D, 1.0_QP, 1.0_QP)
        case (ST_INT)
          p_rnode%p_Iinteger3D = p_rnode%p_Iinteger3D + p_rnodeTmp%p_Iinteger3D
        case (ST_INT8)
          p_rnode%p_Iint8_3D = p_rnode%p_Iint8_3D + p_rnodeTmp%p_Iint8_3D
        case (ST_INT16)
          p_rnode%p_Iint16_3D = p_rnode%p_Iint16_3D + p_rnodeTmp%p_Iint16_3D
        case (ST_INT32)
          p_rnode%p_Iint32_3D = p_rnode%p_Iint32_3D + p_rnodeTmp%p_Iint32_3D
        case (ST_INT64)
          p_rnode%p_Iint64_3D = p_rnode%p_Iint64_3D + p_rnodeTmp%p_Iint64_3D
        case (ST_LOGICAL)
          p_rnode%p_Blogical3D = p_rnode%p_Blogical3D .or. p_rnodeTmp%p_Blogical3D
        case default
          call output_line ('Unsupported data typee!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
          call sys_halt()
        end select
        
      case default
        call output_line ('Unsupported dimension!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
        call sys_halt()
      end select
      
      ! Manually remove memory address from handle ihandleTmp
      p_rnodeTmp%cdeviceMemPtr = C_NULL_PTR

      ! Release temporal handle ihandleTmp
      call storage_free(ihandleTmp)

    case default
      call output_line ('Unsupported synchronisation!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
      call sys_halt()
    end select

#else

    call output_line ('Application must be compiled with coprocessor support enabled!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_syncMemoryHostDevice')
    call sys_halt()

#endif

  contains

    ! Here, the working routines follow

    ! In-place transpose the memory order of a single 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Single2D(Fsingle2D,n1,n2)
      real(SP), dimension(*), intent(inout) :: Fsingle2D
      integer, intent(in) :: n1,n2

      real(SP), dimension(:), allocatable :: Fsingle2D_t

      allocate(Fsingle2D_t(n1*n2))
      Fsingle2D_t(1:n1*n2) = Fsingle2D(1:n1*n2)
      call transposeOutOfPlace_Single2D(Fsingle2D_t,Fsingle2D,n1,n2)
      deallocate(Fsingle2D_t)
            
    end subroutine transposeInPlace_Single2D

    ! In-place transpose the memory order of a single 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Single3D(Fsingle3D,n1,n2,n3)
      real(SP), dimension(*), intent(inout) :: Fsingle3D
      integer, intent(in) :: n1,n2,n3

      real(SP), dimension(:), allocatable :: Fsingle3D_t

      allocate(Fsingle3D_t(n1*n2*n3))
      Fsingle3D_t(1:n1*n2*n3) = Fsingle3D(1:n1*n2*n3)
      call transposeOutOfPlace_Single3D(Fsingle3D_t,Fsingle3D,n1,n2,n3)
      deallocate(Fsingle3D_t)
            
    end subroutine transposeInPlace_Single3D

    ! Out-of-place transpose the memory order of a single 2D array
    ! stored in host memory into a single 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Single2D(Fsingle2D,Fsingle2D_t,n1,n2)
      real(SP), intent(in), dimension(*) :: Fsingle2D
      real(SP), intent(out), dimension(*) :: Fsingle2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Fsingle2D_t(n2*(i-1)+j) = Fsingle2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Single2D

    ! Out-of-place transpose the memory order of a single 3D array
    ! stored in host memory into a single 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Single3D(Fsingle3D,Fsingle3D_t,n1,n2,n3)
      real(SP), intent(in), dimension(*) :: Fsingle3D
      real(SP), intent(out), dimension(*) :: Fsingle3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Fsingle3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Fsingle3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Single3D
    
    ! In-place transpose the memory order of a double 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Double2D(Ddouble2D,n1,n2)
      real(DP), dimension(*), intent(inout) :: Ddouble2D
      integer, intent(in) :: n1,n2

      real(DP), dimension(:), allocatable :: Ddouble2D_t

      allocate(Ddouble2D_t(n1*n2))
      Ddouble2D_t(1:n1*n2) = Ddouble2D(1:n1*n2)
      call transposeOutOfPlace_Double2D(Ddouble2D_t,Ddouble2D,n1,n2)
      deallocate(Ddouble2D_t)
            
    end subroutine transposeInPlace_Double2D

    ! In-place transpose the memory order of a double 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Double3D(Ddouble3D,n1,n2,n3)
      real(DP), dimension(*), intent(inout) :: Ddouble3D
      integer, intent(in) :: n1,n2,n3

      real(DP), dimension(:), allocatable :: Ddouble3D_t

      allocate(Ddouble3D_t(n1*n2*n3))
      Ddouble3D_t(1:n1*n2*n3) = Ddouble3D(1:n1*n2*n3)
      call transposeOutOfPlace_Double3D(Ddouble3D_t,Ddouble3D,n1,n2,n3)
      deallocate(Ddouble3D_t)
            
    end subroutine transposeInPlace_Double3D

    ! Out-of-place transpose the memory order of a double 2D array
    ! stored in host memory into a double 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Double2D(Ddouble2D,Ddouble2D_t,n1,n2)
      real(DP), intent(in), dimension(*) :: Ddouble2D
      real(DP), intent(out), dimension(*) :: Ddouble2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Ddouble2D_t(n2*(i-1)+j) = Ddouble2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Double2D

    ! Out-of-place transpose the memory order of a double 3D array
    ! stored in host memory into a double 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Double3D(Ddouble3D,Ddouble3D_t,n1,n2,n3)
      real(DP), intent(in), dimension(*) :: Ddouble3D
      real(DP), intent(out), dimension(*) :: Ddouble3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Ddouble3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Ddouble3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Double3D

    ! In-place transpose the memory order of a quad 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Quad2D(Qquad2D,n1,n2)
      real(QP), dimension(*), intent(inout) :: Qquad2D
      integer, intent(in) :: n1,n2

      real(QP), dimension(:), allocatable :: Qquad2D_t

      allocate(Qquad2D_t(n1*n2))
      Qquad2D_t(1:n1*n2) = Qquad2D(1:n1*n2)
      call transposeOutOfPlace_Quad2D(Qquad2D_t,Qquad2D,n1,n2)
      deallocate(Qquad2D_t)
            
    end subroutine transposeInPlace_Quad2D

    ! In-place transpose the memory order of a quad 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Quad3D(Qquad3D,n1,n2,n3)
      real(QP), dimension(*), intent(inout) :: Qquad3D
      integer, intent(in) :: n1,n2,n3

      real(QP), dimension(:), allocatable :: Qquad3D_t

      allocate(Qquad3D_t(n1*n2*n3))
      Qquad3D_t(1:n1*n2*n3) = Qquad3D(1:n1*n2*n3)
      call transposeOutOfPlace_Quad3D(Qquad3D_t,Qquad3D,n1,n2,n3)
      deallocate(Qquad3D_t)
            
    end subroutine transposeInPlace_Quad3D

    ! Out-of-place transpose the memory order of a quad 2D array
    ! stored in host memory into a quad 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Quad2D(Qquad2D,Qquad2D_t,n1,n2)
      real(QP), intent(in), dimension(*) :: Qquad2D
      real(QP), intent(out), dimension(*) :: Qquad2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Qquad2D_t(n2*(i-1)+j) = Qquad2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Quad2D

    ! Out-of-place transpose the memory order of a quad 3D array
    ! stored in host memory into a quad 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Quad3D(Qquad3D,Qquad3D_t,n1,n2,n3)
      real(QP), intent(in), dimension(*) :: Qquad3D
      real(QP), intent(out), dimension(*) :: Qquad3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Qquad3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Qquad3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Quad3D

    ! In-place transpose the memory order of an integer 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Integer2D(Iinteger2D,n1,n2)
      integer, dimension(*), intent(inout) :: Iinteger2D
      integer, intent(in) :: n1,n2

      integer, dimension(:), allocatable :: Iinteger2D_t

      allocate(Iinteger2D_t(n1*n2))
      Iinteger2D_t(1:n1*n2) = Iinteger2D(1:n1*n2)
      call transposeOutOfPlace_Integer2D(Iinteger2D_t,Iinteger2D,n1,n2)
      deallocate(Iinteger2D_t)
            
    end subroutine transposeInPlace_Integer2D

    ! In-place transpose the memory order of a integer 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Integer3D(Iinteger3D,n1,n2,n3)
      integer, dimension(*), intent(inout) :: Iinteger3D
      integer, intent(in) :: n1,n2,n3

      integer, dimension(:), allocatable :: Iinteger3D_t
      
      allocate(Iinteger3D_t(n1*n2*n3))
      Iinteger3D_t(1:n1*n2*n3) = Iinteger3D(1:n1*n2*n3)
      call transposeOutOfPlace_Integer3D(Iinteger3D_t,Iinteger3D,n1,n2,n3)
      deallocate(Iinteger3D_t)
            
    end subroutine transposeInPlace_Integer3D

    ! Out-of-place transpose the memory order of a integer 2D array
    ! stored in host memory into a integer 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Integer2D(Iinteger2D,Iinteger2D_t,n1,n2)
      integer, intent(in), dimension(*) :: Iinteger2D
      integer, intent(out), dimension(*) :: Iinteger2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Iinteger2D_t(n2*(i-1)+j) = Iinteger2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Integer2D
    
    ! Out-of-place transpose the memory order of a integer 3D array
    ! stored in host memory into a integer 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Integer3D(Iinteger3D,Iinteger3D_t,n1,n2,n3)
      integer, intent(in), dimension(*) :: Iinteger3D
      integer, intent(out), dimension(*) :: Iinteger3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Iinteger3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Iinteger3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Integer3D

    ! In-place transpose the memory order of an int8 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Int8_2D(Iint8_2D,n1,n2)
      integer(I8), dimension(*), intent(inout) :: Iint8_2D
      integer, intent(in) :: n1,n2

      integer(I8), dimension(:), allocatable :: Iint8_2D_t

      allocate(Iint8_2D_t(n1*n2))
      Iint8_2D_t(1:n1*n2) = Iint8_2D(1:n1*n2)
      call transposeOutOfPlace_Int8_2D(Iint8_2D_t,Iint8_2D,n1,n2)
      deallocate(Iint8_2D_t)
            
    end subroutine transposeInPlace_Int8_2D
    
    ! In-place transpose the memory order of a int8 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Int8_3D(Iint8_3D,n1,n2,n3)
      integer(I8), dimension(*), intent(inout) :: Iint8_3D
      integer, intent(in) :: n1,n2,n3

      integer(I8), dimension(:), allocatable :: Iint8_3D_t

      allocate(Iint8_3D_t(n1*n2*n3))
      Iint8_3D_t(1:n1*n2*n3) = Iint8_3D(1:n1*n2*n3)
      call transposeOutOfPlace_Int8_3D(Iint8_3D_t,Iint8_3D,n1,n2,n3)
      deallocate(Iint8_3D_t)
            
    end subroutine transposeInPlace_Int8_3D

    ! Out-of-place transpose the memory order of a int8 2D array
    ! stored in host memory into a int8 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Int8_2D(Iint8_2D,Iint8_2D_t,n1,n2)
      integer(I8), intent(in), dimension(*) :: Iint8_2D
      integer(I8), intent(out), dimension(*) :: Iint8_2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Iint8_2D_t(n2*(i-1)+j) = Iint8_2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Int8_2D

    ! Out-of-place transpose the memory order of a int8 3D array
    ! stored in host memory into a int8 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Int8_3D(Iint8_3D,Iint8_3D_t,n1,n2,n3)
      integer(I8), intent(in), dimension(*) :: Iint8_3D
      integer(I8), intent(out), dimension(*) :: Iint8_3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Iint8_3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Iint8_3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Int8_3D

    ! In-place transpose the memory order of an int16 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Int16_2D(Iint16_2D,n1,n2)
      integer(I16), dimension(*), intent(inout) :: Iint16_2D
      integer, intent(in) :: n1,n2

      integer(I16), dimension(:), allocatable :: Iint16_2D_t

      allocate(Iint16_2D_t(n1*n2))
      Iint16_2D_t(1:n1*n2) = Iint16_2D(1:n1*n2)
      call transposeOutOfPlace_Int16_2D(Iint16_2D_t,Iint16_2D,n1,n2)
      deallocate(Iint16_2D_t)
            
    end subroutine transposeInPlace_Int16_2D

    ! In-place transpose the memory order of a int16 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Int16_3D(Iint16_3D,n1,n2,n3)
      integer(I16), dimension(*), intent(inout) :: Iint16_3D
      integer, intent(in) :: n1,n2,n3

      integer(I16), dimension(:), allocatable :: Iint16_3D_t

      allocate(Iint16_3D_t(n1*n2*n3))
      Iint16_3D_t(1:n1*n2*n3) = Iint16_3D(1:n1*n2*n3)
      call transposeOutOfPlace_Int16_3D(Iint16_3D_t,Iint16_3D,n1,n2,n3)
      deallocate(Iint16_3D_t)
            
    end subroutine transposeInPlace_Int16_3D

    ! Out-of-place transpose the memory order of a int16 2D array
    ! stored in host memory into a int16 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Int16_2D(Iint16_2D,Iint16_2D_t,n1,n2)
      integer(I16), intent(in), dimension(*) :: Iint16_2D
      integer(I16), intent(out), dimension(*) :: Iint16_2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Iint16_2D_t(n2*(i-1)+j) = Iint16_2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Int16_2D

    ! Out-of-place transpose the memory order of a int16 3D array
    ! stored in host memory into a int16 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Int16_3D(Iint16_3D,Iint16_3D_t,n1,n2,n3)
      integer(I16), intent(in), dimension(*) :: Iint16_3D
      integer(I16), intent(out), dimension(*) :: Iint16_3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Iint16_3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Iint16_3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Int16_3D

    ! In-place transpose the memory order of an int32 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Int32_2D(Iint32_2D,n1,n2)
      integer(I32), dimension(*), intent(inout) :: Iint32_2D
      integer, intent(in) :: n1,n2

      integer(I32), dimension(:), allocatable :: Iint32_2D_t

      allocate(Iint32_2D_t(n1*n2))
      Iint32_2D_t(1:n1*n2) = Iint32_2D(1:n1*n2)
      call transposeOutOfPlace_Int32_2D(Iint32_2D_t,Iint32_2D,n1,n2)
      deallocate(Iint32_2D_t)
            
    end subroutine transposeInPlace_Int32_2D

    ! In-place transpose the memory order of a int32 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Int32_3D(Iint32_3D,n1,n2,n3)
      integer(I32), dimension(*), intent(inout) :: Iint32_3D
      integer, intent(in) :: n1,n2,n3

      integer(I32), dimension(:), allocatable :: Iint32_3D_t

      allocate(Iint32_3D_t(n1*n2*n3))
      Iint32_3D_t(1:n1*n2*n3) = Iint32_3D(1:n1*n2*n3)
      call transposeOutOfPlace_Int32_3D(Iint32_3D_t,Iint32_3D,n1,n2,n3)
      deallocate(Iint32_3D_t)
            
    end subroutine transposeInPlace_Int32_3D

    ! Out-of-place transpose the memory order of a int32 2D array
    ! stored in host memory into a int32 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Int32_2D(Iint32_2D,Iint32_2D_t,n1,n2)
      integer(I32), intent(in), dimension(*) :: Iint32_2D
      integer(I32), intent(out), dimension(*) :: Iint32_2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Iint32_2D_t(n2*(i-1)+j) = Iint32_2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Int32_2D

    ! Out-of-place transpose the memory order of a int32 3D array
    ! stored in host memory into a int32 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Int32_3D(Iint32_3D,Iint32_3D_t,n1,n2,n3)
      integer(I32), intent(in), dimension(*) :: Iint32_3D
      integer(I32), intent(out), dimension(*) :: Iint32_3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Iint32_3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Iint32_3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Int32_3D

    ! In-place transpose the memory order of an int64 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Int64_2D(Iint64_2D,n1,n2)
      integer(I64), dimension(*), intent(inout) :: Iint64_2D
      integer, intent(in) :: n1,n2

      integer(I64), dimension(:), allocatable :: Iint64_2D_t

      allocate(Iint64_2D_t(n1*n2))
      Iint64_2D_t(1:n1*n2) = Iint64_2D(1:n1*n2)
      call transposeOutOfPlace_Int64_2D(Iint64_2D_t,Iint64_2D,n1,n2)
      deallocate(Iint64_2D_t)
            
    end subroutine transposeInPlace_Int64_2D

    ! In-place transpose the memory order of a int64 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Int64_3D(Iint64_3D,n1,n2,n3)
      integer(I64), dimension(*), intent(inout) :: Iint64_3D
      integer, intent(in) :: n1,n2,n3

      integer(I64), dimension(:), allocatable :: Iint64_3D_t

      allocate(Iint64_3D_t(n1*n2*n3))
      Iint64_3D_t(1:n1*n2*n3) = Iint64_3D(1:n1*n2*n3)
      call transposeOutOfPlace_Int64_3D(Iint64_3D_t,Iint64_3D,n1,n2,n3)
      deallocate(Iint64_3D_t)
            
    end subroutine transposeInPlace_Int64_3D

    ! Out-of-place transpose the memory order of a int64 2D array
    ! stored in host memory into a int64 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Int64_2D(Iint64_2D,Iint64_2D_t,n1,n2)
      integer(I64), intent(in), dimension(*) :: Iint64_2D
      integer(I64), intent(out), dimension(*) :: Iint64_2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Iint64_2D_t(n2*(i-1)+j) = Iint64_2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Int64_2D
    
    ! Out-of-place transpose the memory order of a int64 3D array
    ! stored in host memory into a int64 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Int64_3D(Iint64_3D,Iint64_3D_t,n1,n2,n3)
      integer(I64), intent(in), dimension(*) :: Iint64_3D
      integer(I64), intent(out), dimension(*) :: Iint64_3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Iint64_3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Iint64_3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Int64_3D

    ! In-place transpose the memory order of a logical 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Logical2D(Blogical2D,n1,n2)
      logical, dimension(*), intent(inout) :: Blogical2D
      integer, intent(in) :: n1,n2

      logical, dimension(:), allocatable :: Blogical2D_t

      allocate(Blogical2D_t(n1*n2))
      Blogical2D_t(1:n1*n2) = Blogical2D(1:n1*n2)
      call transposeOutOfPlace_Logical2D(Blogical2D_t,Blogical2D,n1,n2)
      deallocate(Blogical2D_t)
            
    end subroutine transposeInPlace_Logical2D

    ! In-place transpose the memory order of a logical 3D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Logical3D(Blogical3D,n1,n2,n3)
      logical, dimension(*), intent(inout) :: Blogical3D
      integer, intent(in) :: n1,n2,n3

      logical, dimension(:), allocatable :: Blogical3D_t

      allocate(Blogical3D_t(n1*n2*n3))
      Blogical3D_t(1:n1*n2*n3) = Blogical3D(1:n1*n2*n3)
      call transposeOutOfPlace_Logical3D(Blogical3D_t,Blogical3D,n1,n2,n3)
      deallocate(Blogical3D_t)
            
    end subroutine transposeInPlace_Logical3D

    ! Out-of-place transpose the memory order of a logical 2D array
    ! stored in host memory into a logical 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Logical2D(Blogical2D,Blogical2D_t,n1,n2)
      logical, intent(in), dimension(*) :: Blogical2D
      logical, intent(out), dimension(*) :: Blogical2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Blogical2D_t(n2*(i-1)+j) = Blogical2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Logical2D

    ! Out-of-place transpose the memory order of a logical 3D array
    ! stored in host memory into a logical 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Logical3D(Blogical3D,Blogical3D_t,n1,n2,n3)
      logical, intent(in), dimension(*) :: Blogical3D
      logical, intent(out), dimension(*) :: Blogical3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Blogical3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Blogical3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Logical3D
    
    ! In-place transpose the memory order of a character 2D array stored
    ! in host memory. In the current implementation, this routine
    ! works out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Char2D(Schar2D,n1,n2)
      character, dimension(*), intent(inout) :: Schar2D
      integer, intent(in) :: n1,n2

      character, dimension(:), allocatable :: Schar2D_t

      allocate(Schar2D_t(n1*n2))
      Schar2D_t(1:n1*n2) = Schar2D(1:n1*n2)
      call transposeOutOfPlace_Char2D(Schar2D_t,Schar2D,n1,n2)
      deallocate(Schar2D_t)
            
    end subroutine transposeInPlace_Char2D

    ! In-place transpose the memory order of a char 3D array stored in
    ! host memory. In the current implementation, this routine works
    ! out-of-place but future implementations will be in-place.
    subroutine transposeInPlace_Char3D(Schar3D,n1,n2,n3)
      character, dimension(*), intent(inout) :: Schar3D
      integer, intent(in) :: n1,n2,n3

      character, dimension(:), allocatable :: Schar3D_t

      allocate(Schar3D_t(n1*n2*n3))
      Schar3D_t(1:n1*n2*n3) = Schar3D(1:n1*n2*n3)
      call transposeOutOfPlace_Char3D(Schar3D_t,Schar3D,n1,n2,n3)
      deallocate(Schar3D_t)
            
    end subroutine transposeInPlace_Char3D
    
    ! Out-of-place transpose the memory order of a char 2D array
    ! stored in host memory into a char 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Char2D(Schar2D,Schar2D_t,n1,n2)
      character, intent(in), dimension(*) :: Schar2D
      character, intent(out), dimension(*) :: Schar2D_t
      integer, intent(in) :: n1,n2

      integer :: i,j

      do j=1,n2
        do i=1,n1
          Schar2D_t(n2*(i-1)+j) = Schar2D(n1*(j-1)+i)
        end do
      end do
    end subroutine transposeOutOfPlace_Char2D

    ! Out-of-place transpose the memory order of a char 3D array
    ! stored in host memory into a char 1D array stored in device
    ! memory
    subroutine transposeOutOfPlace_Char3D(Schar3D,Schar3D_t,n1,n2,n3)
      character, intent(in), dimension(*) :: Schar3D
      character, intent(out), dimension(*) :: Schar3D_t
      integer, intent(in) :: n1,n2,n3

      integer :: i,j,k

      do k=1,n3
        do j=1,n2
          do i=1,n1
            Schar3D_t(n2*n3*(i-1)+n3*(j-1)+k) = Schar3D(n1*n2*(k-1)+n1*(j-1)+i)
          end do
        end do
      end do
    end subroutine transposeOutOfPlace_Char3D

  end subroutine storage_syncMemoryHostDevice

!************************************************************************

!<subroutine>

  subroutine storage_clearMemoryOnDevice (ihandle, rheap)

!<description>
  ! This routine clears a memory block in the device memory
  ! associated with the given handle of the heap. The handle must
  ! already be associated with some memory block.
!</description>

!<input>
  ! Handle of the memory block
  integer, intent(in) :: ihandle
!</input>

!<inputoutput>
  ! OPTIONAL: local heap structure.
  ! If not given, the global heap is used.
  type(t_storageBlock), intent(inout), target, optional :: rheap
!</inputoutput>

!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

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
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_clearMemoryOnDevice')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)

    ! Check if memory address is already in use
    if (.not.storage_isAssociated(p_rnode%cdeviceMemPtr))&
        call coproc_newMemoryOnDevice(p_rnode%cdeviceMemPtr, p_rnode%imemBytes)

    ! Fill memory with zeros
    call coproc_clearMemoryOnDevice(p_rnode%cdeviceMemPtr, p_rnode%imemBytes)

#else

    call output_line ('Application must be compiled with coprocessor support enabled!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_clearMemoryOnDevice')
    call sys_halt()

#endif

  end subroutine storage_clearMemoryOnDevice

!************************************************************************

!<function>

  function storage_getMemPtrOnDevice (ihandle, rheap) result(cdeviceMemPtr)

!<description>
  ! Returns the memory address of the memory block associated with
  ! handle ihandle in device memory.
!</description>

!<input>
  ! Handle of the memory block
  integer, intent(in) :: ihandle

  ! OPTIONAL: local heap structure to initialise. If not given, the
  ! global heap is used.
  type(t_storageBlock), intent(in), target, optional :: rheap
!</input>

!<result>
  ! Memory address of the memory block in device memory
  type(C_PTR) :: cdeviceMemPtr
!</result>

!</function>

#ifdef ENABLE_COPROCESSOR_SUPPORT

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
                        OU_CLASS_ERROR,OU_MODE_STD,'storage_getMemPtrOnDevice')
      call sys_halt()
    end if

    ! Where is the descriptor of the handle?
    p_rnode => p_rheap%p_Rdescriptors(ihandle)
    
    ! Return memory address
    cdeviceMemPtr = p_rnode%cdeviceMemPtr

#else

    call output_line ('Application must be compiled with coprocessor support enabled!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'storage_getMemPtrOnDevice')
    call sys_halt()

    ! Assign zero in this case; this code is not executed at runtime due to the sys_halt() call above,
    ! however, this supresses the following warning from the IFC v12 compiler:
    ! warning #6178: The return value of this FUNCTION has not been defined.   [P_MEMADDRESS]
    cdeviceMemPtr = C_NULL_PTR

#endif

  end function storage_getMemPtrOnDevice

!************************************************************************

!<function>

  function storage_isAssociated(cmemPtr)

!<description>
    ! This function returns .TRUE. if the memory pointer is associated
!</description>

!<input>
    ! Memory pointer
    type(C_PTR), intent(in) :: cmemPtr
!</input>

!<result>
    logical :: storage_isAssociated
!</result>
!</function>

#ifdef HAS_ISO_C_BINDING
    storage_isAssociated = c_associated(cmemPtr)
#else
    storage_isAssociated = (int(cmemPtr%imemAddress,I64) .gt. 0_I64)
#endif
  end function storage_isAssociated

end module storage
