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
  public :: stack_pushback
  public :: stack_backInt
  public :: stack_backSngl
  public :: stack_backDble
  public :: stack_popbackInt
  public :: stack_popbackSngl
  public :: stack_popbackDble

  interface stack_pushback
    module procedure stack_pushbackInt
    module procedure stack_pushbackSngl
    module procedure stack_pushbackDble
  end interface

!<types>

!<typeblock>

  ! Type block for holding a dynamic stack
  type t_stack
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

  subroutine stack_create(rstack,isize,cdataType)

!<description>
    ! Creates a stack with initial memory
!</description>

!<input>
    ! Initial stack size
    integer, intent(IN) :: isize

    ! Stack data type
    integer, intent(IN) :: cdataType
!</input>

!<inputoutput>
    ! Stack
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>
!</subroutine>
    
    select case (cdataType)
    case (ST_INT)
      call storage_new("stack_create","h_StackData",isize,ST_INT,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT)
      
    case (ST_SINGLE)
      call storage_new("stack_create","h_StackData",isize,ST_SINGLE,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT)

    case (ST_DOUBLE)
      call storage_new("stack_create","h_StackData",isize,ST_DOUBLE,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT)
      
    case DEFAULT
      call output_line('Invalid data type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_create')
      call sys_halt()
    end select

    rstack%istackSize=isize
    rstack%idataType=cdataType
  end subroutine stack_create
  
  !************************************************************************

!<subroutine>

  subroutine stack_release(rstack)

!<description>
    ! Release a stack
!</description>

!<inputoutput>
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>
!</subroutine>
    
    if (rstack%h_StackData /= ST_NOHANDLE) call storage_free(rstack%h_StackData)
    
    rstack%istackSize=0
    rstack%istackPosition=0
  end subroutine stack_release
  
  !************************************************************************

!<subroutine>

  subroutine stack_clear(rstack)

!<description>
    ! Clear stack, i.e., reset stack pointer to zero
!</description>

!<inputoutput>
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>
!</subroutine>
    
    rstack%istackPosition=0       
  end subroutine stack_clear

  !************************************************************************

!<function>

  function stack_isempty(rstack) result(bisempty)

!<description>
    ! Returns TRUE if stack is empty
!</description>

!<input>
    type(t_stack), intent(IN) :: rstack
!</input>

!<result>
    logical :: bisempty
!</result>
!</function>

    bisempty=(rstack%istackPosition==0)
  end function stack_isempty

  !************************************************************************

!<function>

  function stack_size(rstack) result(isize)

!<description>
    ! Returns the stack size
!</description>

!<input>
    type(t_stack), intent(IN) :: rstack
!</input>

!<result>
    integer :: isize
!</result>
!</function>

    isize=rstack%istackPosition
  end function stack_size

  !************************************************************************

!<subroutine>
  
  subroutine stack_pushbackInt(rstack,idata)

!<description>
    ! Add an integer value to the top of the stack
!</description>
    
!<input>
    integer, intent(IN) :: idata
!</input>

!<inputoutput>
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: StackData
    
    if (rstack%h_StackData == ST_NOHANDLE) then
      call output_line('Invalid data type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_pushbackInt')
      call sys_halt()
    end if
    
    ! Double storage for stack if required
    if (rstack%istackSize == rstack%istackPosition) then
      call storage_realloc("stack_pushbackInt",2*rstack%istackSize,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT,.true.)
      rstack%istackSize=2*rstack%istackSize
    end if
    
    ! Push to stack
    rstack%istackPosition=rstack%istackPosition+1
    call storage_getbase_int(rstack%h_StackData,StackData)
    StackData(rstack%istackPosition)=idata
  end subroutine stack_pushbackInt

  !************************************************************************

!<subroutine>
  
  subroutine stack_pushbackSngl(rstack,sdata)

!<description>
    ! Add a single value to the top of the stack
!</description>
    
!<input>
    real(SP), intent(IN) :: sdata
!</input>

!<inputoutput>
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    real(SP), dimension(:), pointer :: StackData

    if (rstack%h_StackData == ST_NOHANDLE) then
      call output_line('Invalid data type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_pushbackSngl')
      call sys_halt()
    end if

    ! Double storage for stack if required
    if (rstack%istackSize == rstack%istackPosition) then
      call storage_realloc("stack_pushbackSngl",2*rstack%istackSize,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT,.true.)
      rstack%istackSize=2*rstack%istackSize
    end if

    ! Push to stack
    rstack%istackPosition=rstack%istackPosition+1
    call storage_getbase_single(rstack%h_StackData,StackData)
    StackData(rstack%istackPosition)=sdata
  end subroutine stack_pushbackSngl

  !************************************************************************

!<subroutine>
  
  subroutine stack_pushbackDble(rstack,ddata)

!<description>
    ! Add a double value to the top of the stack
!</description>
    
!<input>
    real(DP), intent(IN) :: ddata
!</input>

!<inputoutput>
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: StackData

    if (rstack%h_StackData == ST_NOHANDLE) then
      call output_line('Invalid data type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_pushbackDble')
      call sys_halt()
    end if

    ! Double storage for stack if required
    if (rstack%istackSize == rstack%istackPosition) then
      call storage_realloc("stack_pushbackDble",2*rstack%istackSize,&
          rstack%h_StackData,ST_NEWBLOCK_NOINIT,.true.)
      rstack%istackSize=2*rstack%istackSize
    end if
    
    ! Push to stack
    rstack%istackPosition=rstack%istackPosition+1
    call storage_getbase_double(rstack%h_StackData,StackData)
    StackData(rstack%istackPosition)=ddata
  end subroutine stack_pushbackDble

  !************************************************************************

!<function>

  function stack_backInt(rstack) result(idata)

!<description>
    ! Return integer value from top of the stack
!</description>

!<input>
    type(t_stack), intent(IN) :: rstack
!</input>

!<result>
    integer :: idata
!</result>
!</function>

    integer, dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_int(rstack%h_StackData,StackData)
      idata=StackData(rstack%istackPosition)
    else
      call output_line('Stack empty!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_backInt')
      call sys_halt()
    end if
    
  end function stack_backInt

  !************************************************************************

!<function>

  function stack_backSngl(rstack) result(sdata)

!<description>
    ! Return single value from top of the stack
!</description>

!<input>
    type(t_stack), intent(IN) :: rstack
!</input>

!<result>
    real(SP) :: sdata
!</result>
!</function>

    real(SP), dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_single(rstack%h_StackData,StackData)
      sdata=StackData(rstack%istackPosition)
    else
      call output_line('Stack empty!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_backSngl')
      call sys_halt()
    end if
    
  end function stack_backSngl

  !************************************************************************

!<function>

  function stack_backDble(rstack) result(ddata)

!<description>
    ! Return a single value from top of the stack
!</description>

!<input>
    type(t_stack), intent(IN) :: rstack
!</input>

!<result>
    real(SP) :: ddata
!</result>
!</function>

    real(DP), dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_double(rstack%h_StackData,StackData)
      ddata=StackData(rstack%istackPosition)
    else
      call output_line('Stack empty!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_backDble')
      call sys_halt()
    end if
    
  end function stack_backDble

  !************************************************************************

!<function>

  function stack_popbackInt(rstack) result(idata)

!<description>
    ! Remove an integer value from top of the stack
!</description>

!<inputoutput>
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>

!<result>
    integer :: idata
!</result>
!</function>

    integer, dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_int(rstack%h_StackData,StackData)
      idata=StackData(rstack%istackPosition)
      rstack%istackPosition=rstack%istackPosition-1
    else
      call output_line('Stack empty!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_popbackInt')
      call sys_halt()
    end if
    
  end function stack_popbackInt

  !************************************************************************

!<function>

  function stack_popbackSngl(rstack) result(sdata)

!<description>
    ! Remove a single value from top of the stack
!</description>

!<inputoutput>
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>

!<result>
    real(SP) :: sdata
!</result>
!</function>

    real(SP), dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_single(rstack%h_StackData,StackData)
      sdata=StackData(rstack%istackPosition)
      rstack%istackPosition=rstack%istackPosition-1
    else
      call output_line('Stack empty!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_popbackSngl')
      call sys_halt()
    end if
    
  end function stack_popbackSngl

  !************************************************************************

!<function>

  function stack_popbackDble(rstack) result(ddata)

!<description>
    ! Remove a single value from top of the stack
!</description>

!<inputoutput>
    type(t_stack), intent(INOUT) :: rstack
!</inputoutput>

!<result>
    real(SP) :: ddata
!</result>
!</function>

    real(DP), dimension(:), pointer :: StackData
    
    if (.not.stack_isempty(rstack)) then
      call storage_getbase_double(rstack%h_StackData,StackData)
      ddata=StackData(rstack%istackPosition)
      rstack%istackPosition=rstack%istackPosition-1
    else
      call output_line('Stack empty!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stack_popbackDble')
      call sys_halt()
    end if
    
  end function stack_popbackDble
end module stack
