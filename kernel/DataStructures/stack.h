!-*- mode: f90; -*-

#ifndef _STACK_H_
#define _STACK_H_

!##############################################################################
!# ****************************************************************************
!# <name> TEMPLATE(stack,T) </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module header file implements a dynamic stack for T_TYPE data
!#
!# The following routines are available:
!#
!# 1.) stack_create
!#     -> Creates a new empty stack
!#
!# 2.) stack_release
!#     -> Releases an existing stack
!#
!# 3.) stack_clear
!#     -> Clears an existing stack
!#
!# 4.) stack_isempty
!#     -> Tests if an existing stack is empty
!#
!# 5.) stack_size
!#     -> Returns the size of the stack
!#
!# 6.) stack_top
!#     -> Returns the entry from the top of the stack
!#
!# 7.) stack_push
!#     -> Pushes a new entry on top of the stack
!#
!# 8.) stack_pop
!#     -> Pops the entry from the top of the stack
!#
!# </purpose>
!##############################################################################

#include "template.h"

module TEMPLATE(stack,T)

!$use omp_lib
  use fsystem
  use genoutput
  use storage

#ifndef T_STORAGE
  __external_use__(T_MODULE)
#endif
  
  implicit none

  private
  public :: TEMPLATE(t_stack,T)
  public :: stack_create
  public :: stack_release
  public :: stack_clear
  public :: stack_isempty
  public :: stack_size
  public :: stack_top
  public :: stack_push
  public :: stack_pop

  interface stack_create
    module procedure TEMPLATE(stack_create,T)
  end interface

  interface stack_release
    module procedure TEMPLATE(stack_release,T)
  end interface

  interface stack_clear
    module procedure TEMPLATE(stack_clear,T)
  end interface

  interface stack_isempty
    module procedure TEMPLATE(stack_isempty,T)
  end interface

  interface stack_size
    module procedure TEMPLATE(stack_size,T)
  end interface

  interface stack_top
    module procedure TEMPLATE(stack_top,T)
  end interface

  interface stack_push
    module procedure TEMPLATE(stack_push,T)
  end interface
  
  interface stack_pop
    module procedure TEMPLATE(stack_pop,T)
  end interface

!<types>

!<typeblock>

  ! Type block for holding a dynamic stack
  type TEMPLATE(t_stack,T)
    private

    ! Size of stack
    integer :: istackSize = 0

    ! Position of last stack item
    integer :: istackPosition = 0

#ifdef T_STORAGE
    ! Handle to stack data
    integer :: h_StackData = ST_NOHANDLE
#endif

    ! Pointer to stack data
    TTYPE(T_TYPE), dimension(:), pointer :: p_StackData => null()

  end type

!</typeblock>

!</types>

contains

  !************************************************************************

!<subroutine>

  subroutine TEMPLATE(stack_create,T)(rstack, isize)

!<description>
    ! Creates a stack with prescribed initial memory
!</description>

!<input>
    ! Initial stack size
    integer, intent(in) :: isize
!</input>

!<output>
    ! Stack
    type(TEMPLATE(t_stack,T)), intent(out) :: rstack
!</output>
!</subroutine>

    ! Set stack size
    rstack%istackSize = isize

#ifdef T_STORAGE
    call storage_new('stack_create','h_StackData',rstack%istackSize,&
        T_STORAGE,rstack%h_StackData,ST_NEWBLOCK_NOINIT)
    call storage_getbase(rstack%h_StackData,rstack%p_StackData)
#else
    allocate(rstack%p_StackData(rstack%istackSize))
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine TEMPLATE(stack_release,T)(rstack)

!<description>
    ! Release a stack
!</description>

!<inputoutput>
    ! Stack
    type(TEMPLATE(t_stack,T)), intent(inout) :: rstack
!</inputoutput>
!</subroutine>
    
#ifdef T_STORAGE
    if (rstack%h_StackData .ne. ST_NOHANDLE)&
        call storage_free(rstack%h_StackData)
    nullify(rstack%p_StackData)
#else
    deallocate(rstack%p_StackData)
#endif

    rstack%istackSize     = 0
    rstack%istackPosition = 0

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine TEMPLATE(stack_clear,T)(rstack)

!<description>
    ! Clear stack, i.e., reset stack pointer to zero
!</description>

!<inputoutput>
    ! Stack
    type(TEMPLATE(t_stack,T)), intent(inout) :: rstack
!</inputoutput>
!</subroutine>
    
    rstack%istackPosition = 0

  end subroutine

  !************************************************************************

!<function>

  function TEMPLATE(stack_isEmpty,T)(rstack) result(bisempty)

!<description>
    ! Checks if the stack is empty
!</description>

!<input>
    ! Stack
    type(TEMPLATE(t_stack,T)), intent(in) :: rstack
!</input>

!<result>
    ! Flag: is true if the stack is empty
    logical :: bisempty
!</result>
!</function>

    bisEmpty = (rstack%istackPosition .eq. 0)

  end function

  !************************************************************************

!<function>

  function TEMPLATE(stack_size,T)(rstack) result(isize)

!<description>
    ! Returns the stack size
!</description>

!<input>
    ! Stack
    type(TEMPLATE(t_stack,T)), intent(in) :: rstack
!</input>

!<result>
    ! Size of the stack
    integer :: isize
!</result>
!</function>

    isize = rstack%istackPosition

  end function

  !************************************************************************

!<subroutine>
  
  subroutine TEMPLATE(stack_push,T)(rstack, data)

!<description>
    ! Add a new value to the top of the stack
!</description>
    
!<input>
    ! Data to be pushed onto the stack
    TTYPE(T_TYPE), intent(in) :: data
!</input>

!<inputoutput>
    ! Stack
    type(TEMPLATE(t_stack,T)), intent(inout) :: rstack
!</inputoutput>
!</subroutine>
    
#ifdef T_STORAGE
    if (rstack%h_StackData .eq. ST_NOHANDLE) then
      call output_line('Invalid data type!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_push')
      call sys_halt()
    end if
    
    ! Double storage for stack if required
    if (rstack%istackSize .eq. rstack%istackPosition) then
      rstack%istackSize = 2*rstack%istackSize
      call storage_realloc('stack_push', rstack%istackSize,&
          rstack%h_StackData, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase(rstack%h_StackData, rstack%p_StackData)
    end if
#else
    ! local variable
    type(T_TYPE), dimension(:), pointer :: p_StackData
    
    ! Double storage for stack if required
    if (rstack%istackSize .eq. rstack%istackPosition) then
      allocate(p_StackData(rstack%istackSize))
      p_StackData = rstack%p_StackData
      deallocate(rstack%p_StackData)
      rstack%istackSize = 2*rstack%istackSize
      allocate(rstack%p_StackData(rstack%istackSize))
      rstack%p_StackData(1:size(p_StackData)) = p_StackData
      deallocate(p_StackData)
    end if
#endif
    
    ! Push data to stack
    rstack%istackPosition = rstack%istackPosition+1
    rstack%p_StackData(rstack%istackPosition) = data

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine TEMPLATE(stack_top,T)(rstack, data)

!<description>
    ! Return value from top of the stack
!</description>

!<input>
    ! Stack
    type(TEMPLATE(t_stack,T)), intent(in) :: rstack
!</input>

!<output>
    ! Data on top of the stack
    TTYPE(T_TYPE), intent(out) :: data
!</output>
!</subroutine>

    if (.not.stack_isempty(rstack)) then
      data = rstack%p_StackData(rstack%istackPosition)
    else
      call output_line('Stack empty!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_top')
      call sys_halt()
    end if
    
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine TEMPLATE(stack_pop,T)(rstack, data)

!<description>
    ! Remove a value from the top of the stack
!</description>

!<inputoutput>
    ! Stack
    type(TEMPLATE(t_stack,T)), intent(inout) :: rstack
!</inputoutput>

!<output>
    ! Item removed from the top of the stack
    TTYPE(T_TYPE), intent(out) :: data
!</output>
!</subroutine>
    
    call stack_top(rstack, data)
    rstack%istackPosition = rstack%istackPosition-1

  end subroutine

end module

#endif
