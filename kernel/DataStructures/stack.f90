!##############################################################################
!# ****************************************************************************
!# <name> stack </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a dynamic stack which can be use to store
!# integer, single and double data.
!#
!# </purpose>
!##############################################################################

module stack

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  
  implicit none

  private

  public :: t_stack
  public :: stack_create
  public :: stack_release
  public :: stack_clear
  public :: stack_isempty
  public :: stack_size
  public :: stack_top
  public :: stack_topInt
  public :: stack_topSngl
  public :: stack_topDble
  public :: stack_push
  public :: stack_pushInt
  public :: stack_pushSngl
  public :: stack_pushDble
  public :: stack_pop
  public :: stack_popInt
  public :: stack_popSngl
  public :: stack_popDble

  interface stack_top
    module procedure stack_topInt
    module procedure stack_topSngl
    module procedure stack_topDble
  end interface

  interface stack_push
    module procedure stack_pushInt
    module procedure stack_pushSngl
    module procedure stack_pushDble
  end interface

  interface stack_pop
    module procedure stack_popInt
    module procedure stack_popSngl
    module procedure stack_popDble
  end interface

!<types>

!<typeblock>

  ! Type block for holding a dynamic stack
  type t_stack
    private

    ! Size of stack
    integer :: istackSize = 0

    ! Position of last stack item
    integer :: istackPosition = 0

    ! Data type stored in stack
    integer :: idataType

    ! Handle to stack data
    integer :: h_StackData = ST_NOHANDLE
  end type t_stack

!</typeblock>

!</types>
  
contains

!<subroutine>

  subroutine stack_create(rstack, isize, cdataType)

!<description>
    ! Creates a stack with prescribed initial memory
!</description>

!<input>
    ! Initial stack size
    integer, intent(in) :: isize

    ! Stack data type
    integer, intent(in) :: cdataType
!</input>

!<inputoutput>
    ! Stack
    type(t_stack), intent(inout) :: rstack
!</inputoutput>
!</subroutine>
    
    select case (cdataType)
    case (ST_INT)
      call storage_new('stack_create','h_StackData',isize,ST_INT,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT)
      
    case (ST_SINGLE)
      call storage_new('stack_create','h_StackData',isize,ST_SINGLE,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT)

    case (ST_DOUBLE)
      call storage_new('stack_create','h_StackData',isize,ST_DOUBLE,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT)
      
    case DEFAULT
      call output_line('Invalid data type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_create')
      call sys_halt()
    end select

    rstack%istackSize = isize
    rstack%idataType = cdataType

  end subroutine stack_create
  
  !************************************************************************

!<subroutine>

  subroutine stack_release(rstack)

!<description>
    ! Release a stack
!</description>

!<inputoutput>
    type(t_stack), intent(inout) :: rstack
!</inputoutput>
!</subroutine>
    
    if (rstack%h_StackData .ne. ST_NOHANDLE)&
        call storage_free(rstack%h_StackData)
    
    rstack%istackSize = 0
    rstack%istackPosition = 0

  end subroutine stack_release
  
  !************************************************************************

!<subroutine>

  subroutine stack_clear(rstack)

!<description>
    ! Clear stack, i.e., reset stack pointer to zero
!</description>

!<inputoutput>
    type(t_stack), intent(inout) :: rstack
!</inputoutput>
!</subroutine>
    
    rstack%istackPosition = 0

  end subroutine stack_clear

  !************************************************************************

!<function>

  function stack_isEmpty(rstack) result(bisempty)

!<description>
    ! Returns TRUE if stack is empty
!</description>

!<input>
    type(t_stack), intent(in) :: rstack
!</input>

!<result>
    logical :: bisempty
!</result>
!</function>

    bisEmpty = (rstack%istackPosition .eq. 0)

  end function stack_isEmpty

  !************************************************************************

!<function>

  function stack_size(rstack) result(isize)

!<description>
    ! Returns the stack size
!</description>

!<input>
    type(t_stack), intent(in) :: rstack
!</input>

!<result>
    integer :: isize
!</result>
!</function>

    isize = rstack%istackPosition

  end function stack_size

  !************************************************************************

!<subroutine>
  
  subroutine stack_pushInt(rstack, idata)

!<description>
    ! Add an integer value to the top of the stack
!</description>
    
!<input>
    integer, intent(in) :: idata
!</input>

!<inputoutput>
    type(t_stack), intent(inout) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: StackData
    
    if (rstack%h_StackData .eq. ST_NOHANDLE) then
      call output_line('Invalid data type!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_pushInt')
      call sys_halt()
    end if
    
    ! Double storage for stack if required
    if (rstack%istackSize .eq. rstack%istackPosition) then
      call storage_realloc('stack_pushInt', 2*rstack%istackSize,&
          rstack%h_StackData, ST_NEWBLOCK_NOINIT, .true.)
      rstack%istackSize = 2*rstack%istackSize
    end if
    
    ! Push to stack
    rstack%istackPosition = rstack%istackPosition+1
    call storage_getbase_int(rstack%h_StackData, StackData)
    StackData(rstack%istackPosition) = idata

  end subroutine stack_pushInt

  !************************************************************************

!<subroutine>
  
  subroutine stack_pushSngl(rstack, sdata)

!<description>
    ! Add a single value to the top of the stack
!</description>
    
!<input>
    real(SP), intent(in) :: sdata
!</input>

!<inputoutput>
    type(t_stack), intent(inout) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    real(SP), dimension(:), pointer :: StackData

    if (rstack%h_StackData .eq. ST_NOHANDLE) then
      call output_line('Invalid data type!',&
          OU_CLASS_ERROR, OU_MODE_STD,' stack_pushSngl')
      call sys_halt()
    end if

    ! Double storage for stack if required
    if (rstack%istackSize .eq. rstack%istackPosition) then
      call storage_realloc('stack_pushSngl', 2*rstack%istackSize,&
          rstack%h_StackData, ST_NEWBLOCK_NOINIT, .true.)
      rstack%istackSize = 2*rstack%istackSize
    end if

    ! Push to stack
    rstack%istackPosition = rstack%istackPosition+1
    call storage_getbase_single(rstack%h_StackData, StackData)
    StackData(rstack%istackPosition) = sdata

  end subroutine stack_pushSngl

  !************************************************************************

!<subroutine>
  
  subroutine stack_pushDble(rstack, ddata)

!<description>
    ! Add a double value to the top of the stack
!</description>
    
!<input>
    real(DP), intent(in) :: ddata
!</input>

!<inputoutput>
    type(t_stack), intent(inout) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: StackData

    if (rstack%h_StackData .eq. ST_NOHANDLE) then
      call output_line('Invalid data type!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_pushDble')
      call sys_halt()
    end if

    ! Double storage for stack if required
    if (rstack%istackSize .eq. rstack%istackPosition) then
      call storage_realloc('stack_pushDble', 2*rstack%istackSize,&
          rstack%h_StackData, ST_NEWBLOCK_NOINIT, .true.)
      rstack%istackSize = 2*rstack%istackSize
    end if
    
    ! Push to stack
    rstack%istackPosition = rstack%istackPosition+1
    call storage_getbase_double(rstack%h_StackData, StackData)
    StackData(rstack%istackPosition) = ddata

  end subroutine stack_pushDble

  !************************************************************************

!<subroutine>

  subroutine stack_topInt(rstack, idata)

!<description>
    ! Return integer value from top of the stack
!</description>

!<input>
    type(t_stack), intent(in) :: rstack
!</input>

!<inputoutput>
    integer, intent(inout) :: idata
!</inputoutput>
!</subroutine>

    integer, dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_int(rstack%h_StackData, StackData)
      idata = StackData(rstack%istackPosition)
    else
      call output_line('Stack empty!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_topInt')
      call sys_halt()
    end if
    
  end subroutine stack_topInt

  !************************************************************************

!<subroutine>

  subroutine stack_topSngl(rstack, sdata)

!<description>
    ! Return single value from top of the stack
!</description>

!<input>
    type(t_stack), intent(in) :: rstack
!</input>

!<inputoutput>
    real(SP) :: sdata
!</inputoutput>
!</subroutine>

    real(SP), dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_single(rstack%h_StackData, StackData)
      sdata = StackData(rstack%istackPosition)
    else
      call output_line('Stack empty!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_topSngl')
      call sys_halt()
    end if
    
  end subroutine stack_topSngl

  !************************************************************************

!<subroutine>

  subroutine stack_topDble(rstack, ddata)

!<description>
    ! Return a single value from top of the stack
!</description>

!<input>
    type(t_stack), intent(in) :: rstack
!</input>

!<inputoutput>
    real(DP) :: ddata
!</inputoutput>
!</subroutine>

    real(DP), dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_double(rstack%h_StackData, StackData)
      ddata = StackData(rstack%istackPosition)
    else
      call output_line('Stack empty!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_topDble')
      call sys_halt()
    end if
    
  end subroutine stack_topDble

  !************************************************************************

!<subroutine>

  subroutine stack_popInt(rstack, idata)

!<description>
    ! Remove an integer value from top of the stack
!</description>

!<inputoutput>
    type(t_stack), intent(inout) :: rstack
    integer, intent(inout) :: idata
!</inputoutput>
!</subroutine>

    integer, dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_int(rstack%h_StackData, StackData)
      idata = StackData(rstack%istackPosition)
      rstack%istackPosition = rstack%istackPosition-1
    else
      call output_line('Stack empty!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_popInt')
      call sys_halt()
    end if
    
  end subroutine stack_popInt

  !************************************************************************

!<subroutine>

  subroutine stack_popSngl(rstack, sdata)

!<description>
    ! Remove a single value from top of the stack
!</description>

!<inputoutput>
    type(t_stack), intent(inout) :: rstack
    real(SP), intent(inout) :: sdata
!</inputoutput>
!</subroutine>

    real(SP), dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_single(rstack%h_StackData, StackData)
      sdata = StackData(rstack%istackPosition)
      rstack%istackPosition = rstack%istackPosition-1
    else
      call output_line('Stack empty!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_popSngl')
      call sys_halt()
    end if
    
  end subroutine stack_popSngl

  !************************************************************************

!<subroutine>

  subroutine stack_popDble(rstack, ddata)

!<description>
    ! Remove a single value from top of the stack
!</description>

!<inputoutput>
    type(t_stack), intent(inout) :: rstack
    real(DP), intent(inout) :: ddata
!</inputoutput>
!</subroutine>

    real(DP), dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_double(rstack%h_StackData, StackData)
      ddata = StackData(rstack%istackPosition)
      rstack%istackPosition = rstack%istackPosition-1
    else
      call output_line('Stack empty!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_popDble')
      call sys_halt()
    end if
    
  end subroutine stack_popDble
end module stack
