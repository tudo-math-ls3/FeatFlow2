#ifndef _STACK_H_
#define _STACK_H_

!##############################################################################
!# ****************************************************************************
!# <name> FEAT2_PP_TEMPLATE_T(stack,T) </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module header file implements a LIFO stack.
!#
!# Stacks are a type of container adaptor, specifically designed to
!# operate in a LIFO context (last-in first-out), where elements are
!# inserted and extracted only from the end of the container. Elements
!# are pushed/popped from the "back" of the specific container, which
!# is known as the top of the stack.
!#
!# -----------------------------------------------------------------------------
!# This module is designed with greatest compatibility to the C++
!# Standard Template Library (STL) concerning naming conventions.
!# -----------------------------------------------------------------------------
!#
!# The following routines are available:
!#
!# 1.) stack_create (constructor in STL)
!#     -> Creates a new empty stack
!#
!# 2.) stack_release (destructor in STL)
!#     -> Releases an existing stack
!#
!# 3.) stack_clear (no equivalent in STL)
!#     -> Clears an existing stack
!#
!# 4.) stack_empty (empty in STL)
!#     -> Tests if the stack is empty
!#
!# 5.) stack_size (size in STL)
!#     -> Returns the size of the stack
!#
!# 6.) stack_top (top in STL)
!#     -> Returns the entry from the top of the stack
!#
!# 7.) stack_push (push in STL)
!#     -> Pushes a new entry on top of the stack
!#
!# 8.) stack_pop (pop in STL)
!#     -> Pops the entry from the top of the stack
!#
!# 9.) stack_contains (no equivalent in STL)
!#     -> Checks if the stack contains a given item
!#
!# 10.) stack_cast
!#      -> Casts a stack to a generic object
!#
!# 11.) stack_uncast
!#      -> Casts a generic object to a stack
!#
!#
!# The following operators are available:
!#
!# 1.) "="  assignment operator
!#
!# 2.) "==" Compares two stacks for equality. Two stacks are equal if
!#          they contain the same number of elements and if they are
!#          equal element-by-element.
!#
!# 3.) "/=" Compares two stacks for non-equality. Two stacks are not equal
!#          if they contain a different number of elements or if they are
!#          not equal element-by-element.
!#
!# 4.) "<"  Lexicographical ordering of two stacks.
!#
!# 5.) "<=" Lexicographical ordering of two stacks.
!#
!# 6.) ">"  Lexicographical ordering of two stacks.
!#
!# 7.) ">=" Lexicographical ordering of two stacks.
!#
!# </purpose>
!##############################################################################

#include "../template.h"

  implicit none

  private
  public :: FEAT2_PP_TEMPLATE_T(t_stack,T)
  public :: stack_create
  public :: stack_release
  public :: stack_clear
  public :: stack_empty
  public :: stack_size
  public :: stack_top
  public :: stack_push
  public :: stack_pop
  public :: stack_contains
  public :: stack_cast
  public :: stack_uncast

  public assignment(=)
  public operator(==)
  public operator(/=)
  public operator(<)
  public operator(<=)
  public operator(>)
  public operator(>=)

  interface stack_create
    module procedure FEAT2_PP_TEMPLATE_T(stack_create,T)
  end interface

  interface stack_release
    module procedure FEAT2_PP_TEMPLATE_T(stack_release,T)
  end interface

  interface stack_clear
    module procedure FEAT2_PP_TEMPLATE_T(stack_clear,T)
  end interface

  interface stack_empty
    module procedure FEAT2_PP_TEMPLATE_T(stack_empty,T)
  end interface

  interface stack_size
    module procedure FEAT2_PP_TEMPLATE_T(stack_size,T)
  end interface

  interface stack_top
    module procedure FEAT2_PP_TEMPLATE_T(stack_top,T)
  end interface

  interface stack_push
    module procedure FEAT2_PP_TEMPLATE_T(stack_push,T)
  end interface

  interface stack_pop
    module procedure FEAT2_PP_TEMPLATE_T(stack_pop,T)
  end interface

  interface stack_contains
    module procedure FEAT2_PP_TEMPLATE_T(stack_contains,T)
  end interface

  interface stack_cast
    module procedure FEAT2_PP_TEMPLATE_T(stack_cast,T)
  end interface

  interface stack_uncast
    module procedure FEAT2_PP_TEMPLATE_T(stack_uncast,T)
  end interface

  interface assignment(=)
    module procedure FEAT2_PP_TEMPLATE_T(stack_fassign,T)
  end interface

  interface operator(==)
    module procedure FEAT2_PP_TEMPLATE_T(stack_eq,T)
  end interface

  interface operator(/=)
    module procedure FEAT2_PP_TEMPLATE_T(stack_ne,T)
  end interface

  interface operator(<)
    module procedure FEAT2_PP_TEMPLATE_T(stack_lt,T)
  end interface

  interface operator(<=)
    module procedure FEAT2_PP_TEMPLATE_T(stack_le,T)
  end interface

  interface operator(>)
    module procedure FEAT2_PP_TEMPLATE_T(stack_gt,T)
  end interface

  interface operator(>=)
    module procedure FEAT2_PP_TEMPLATE_T(stack_ge,T)
  end interface

!<types>

!<typeblock>

  ! Type block for holding a dynamic stack
  type FEAT2_PP_TEMPLATE_T(t_stack,T)
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
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), pointer :: p_StackData => null()

  end type

!</typeblock>

!</types>

contains

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(stack_create,T)(rstack, isize)

!<description>
    ! Creates a stack with prescribed initial memory
!</description>

!<input>
    ! Initial stack size
    integer, intent(in) :: isize
!</input>

!<output>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(out) :: rstack
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

  subroutine FEAT2_PP_TEMPLATE_T(stack_release,T)(rstack)

!<description>
    ! Release a stack
!</description>

!<inputoutput>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(inout) :: rstack
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    if (rstack%h_StackData .ne. ST_NOHANDLE)&
        call storage_free(rstack%h_StackData)
    nullify(rstack%p_StackData)
#else
    if (associated(rstack%p_StackData))&
        deallocate(rstack%p_StackData)
#endif

    rstack%istackSize     = 0
    rstack%istackPosition = 0

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(stack_clear,T)(rstack)

!<description>
    ! Clear stack, i.e., reset stack pointer to zero
!</description>

!<inputoutput>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(inout) :: rstack
!</inputoutput>
!</subroutine>

    rstack%istackPosition = 0

  end subroutine

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(stack_Empty,T)(rstack) result(bempty)

!<description>
    ! Checks if the stack is empty
!</description>

!<input>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack
!</input>

!<result>
    ! Flag: is true if the stack is empty
    logical :: bempty
!</result>
!</function>

    bempty = (rstack%istackPosition .eq. 0)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(stack_size,T)(rstack) result(isize)

!<description>
    ! Returns the stack size
!</description>

!<input>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack
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

  subroutine FEAT2_PP_TEMPLATE_T(stack_push,T)(rstack, data)

!<description>
    ! Add a new value to the top of the stack
!</description>

!<input>
    ! Data to be pushed onto the stack
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: data
!</input>

!<inputoutput>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(inout) :: rstack
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

  subroutine FEAT2_PP_TEMPLATE_T(stack_top,T)(rstack, data)

!<description>
    ! Return value from top of the stack
!</description>

!<input>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack
!</input>

!<output>
    ! Data on top of the stack
    FEAT2_PP_TTYPE(T_TYPE), intent(out) :: data
!</output>
!</subroutine>

    if (.not.stack_empty(rstack)) then
      data = rstack%p_StackData(rstack%istackPosition)
    else
      call output_line('Stack empty!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'stack_top')
      call sys_halt()
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(stack_pop,T)(rstack, data)

!<description>
    ! Remove a value from the top of the stack
!</description>

!<inputoutput>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(inout) :: rstack
!</inputoutput>

!<output>
    ! Item removed from the top of the stack
    FEAT2_PP_TTYPE(T_TYPE), intent(out) :: data
!</output>
!</subroutine>

    call stack_top(rstack, data)
    rstack%istackPosition = rstack%istackPosition-1

  end subroutine

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(stack_contains,T)(rstack, data) result(bresult)

!<description>
    ! Check if the stack contains the given data item
!</description>

!<input>
    ! Stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack

    ! Item to be searched for in the stack
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: data
!</inputoutput>

!<result>
    ! True if stack contains the given data item
    logical :: bresult
!</result>
!</function>

    ! local variable
    integer :: i

    bresult = .false.

    do i = rstack%istackPosition, 1, -1
      bresult = (rstack%p_StackData(i) .eq. data)
      if (bresult) return
    end do

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(stack_fassign,T)(rstackDest,rstackSrc)

!<description>
    ! Assigns the content of rstackSrc to rstackDest
!</description>

!<input>
    ! Source stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstackSrc
!</input>

!<output>
    ! Destination stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(out) :: rstackDest
!</output>
!</subroutine>

    ! local variable
    integer :: i

    ! Create empty stack
    call stack_create(rstackDest, rstackSrc%istackSize)

    do i = 1, rstackSrc%istackSize
      rstackDest%p_StackData(i) = rstackSrc%p_StackData(i)
    end do

  end subroutine

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(stack_eq,T)(rstack1,rstack2) result(beq)

!<description>
    ! Compare two stacks for equality
!</description>

!<input>
    ! Stacks
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack1,rstack2
!</input>

!<result>
    ! True if both stacks contain the same number of elements and if
    ! they are equal element-by-element
    logical :: beq
!</result>
!</function>

    ! local variable
    integer :: i

    ! Initialisation
    beq = (rstack1%istackSize == rstack2%istackSize)

    ! Early return?
    if (.not.beq) return

    do i = 1, rstack1%istackSize
      beq = (rstack1%p_StackData(i) == rstack2%p_StackData(i))
      if (.not.beq) return
    end do

    ! If we are here, then both stacks are equal

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(stack_ne,T)(rstack1,rstack2) result(bne)

!<description>
    ! Compare two stacks for non-equality
!</description>

!<input>
    ! Stacks
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack1,rstack2
!</input>

!<result>
    ! True if both stacks contain a different same number of elements
    ! or if they are not equal element-by-element
    logical :: bne
!</result>
!</function>

    ! local variable
    integer :: i

    ! Initialisation
    bne = (rstack1%istackSize /= rstack2%istackSize)

    ! Early return?
    if (bne) return

    do i = 1, rstack1%istackSize
      bne = (rstack1%p_StackData(i) /= rstack2%p_StackData(i))
      if (bne) return
    end do

    ! If we are here, then both stacks are equal

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(stack_lt,T)(rstack1,rstack2) result(blt)

!<description>
    ! Checks lexicographical ordering of two stacks
!</description>

!<input>
    ! Stacks
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack1,rstack2
!</input>

!<result>
    ! True if the lexicographical ordering of stack1 is smaller than
    ! that of stack2
    logical :: blt
!</result>
!</function>

    ! local variable
    integer :: i

    ! Initialisation
    blt = (rstack1%istackSize < rstack2%istackSize)

    ! Early return?
    if (.not.blt) return

    do i = 1, min(rstack1%istackSize,rstack2%istackSize)
      blt = (rstack1%p_StackData(i) < rstack2%p_StackData(i))
      if (.not.blt) return
    end do

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(stack_le,T)(rstack1,rstack2) result(ble)

!<description>
    ! Checks lexicographical ordering of two stacks
!</description>

!<input>
    ! Stacks
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack1,rstack2
!</input>

!<result>
    ! True if the lexicographical ordering of stack1 is smaller than
    ! or equal to that of stack2
    logical :: ble
!</result>
!</function>

    ! local variable
    integer :: i

    ! Initialisation
    ble = (rstack1%istackSize <= rstack2%istackSize)

    ! Early return?
    if (.not.ble) return

    do i = 1, min(rstack1%istackSize,rstack2%istackSize)
      ble = (rstack1%p_StackData(i) <= rstack2%p_StackData(i))
      if (.not.ble) return
    end do

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(stack_gt,T)(rstack1,rstack2) result(bgt)

!<description>
    ! Checks lexicographical ordering of two stacks
!</description>

!<input>
    ! Stacks
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack1,rstack2
!</input>

!<result>
    ! True if the lexicographical ordering of stack1 is greater than
    ! that of stack2
    logical :: bgt
!</result>
!</function>

    ! local variable
    integer :: i

    ! Initialisation
    bgt = (rstack1%istackSize > rstack2%istackSize)

    ! Early return?
    if (.not.bgt) return

    do i = 1, min(rstack1%istackSize,rstack2%istackSize)
      bgt = (rstack1%p_StackData(i) > rstack2%p_StackData(i))
      if (.not.bgt) return
    end do

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(stack_ge,T)(rstack1,rstack2) result(bge)

!<description>
    ! Checks lexicographical ordering of two stacks
!</description>

!<input>
    ! Stacks
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in) :: rstack1,rstack2
!</input>

!<result>
    ! True if the lexicographical ordering of stack1 is greater than
    ! or equal to that of stack2
    logical :: bge
!</result>
!</function>

    ! local variable
    integer :: i

    ! Initialisation
    bge = (rstack1%istackSize >= rstack2%istackSize)

    ! Early return?
    if (.not.bge) return

    do i = 1, min(rstack1%istackSize,rstack2%istackSize)
      bge = (rstack1%p_StackData(i) >= rstack2%p_StackData(i))
      if (.not.bge) return
    end do

  end function
  
  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(stack_cast,T,D)(rstack, rgenericObject)

!<description>
    ! This subroutine casts the given stack to a generic object.
!</description>

!<input>
    ! The stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), intent(in), target :: rstack
!</input>

!<output>
    ! The generic object
    type(t_genericObject), intent(out) :: rgenericObject
!</output>
!</subroutine>

    ! Internal data structure
    type t_void_ptr
      type(FEAT2_PP_TEMPLATE_T(t_stack,T)), pointer :: p_robj => null()
    end type t_void_ptr
    
    ! Internal variables
    type(t_void_ptr) :: rptr

    ! Wrap stack by void pointer structure
    rptr%p_robj => rstack
    
    ! Determine the size of the void pointer structure
    rgenericObject%isize = size(transfer(rptr, rgenericObject%p_cdata))
    
    ! Allocate memory and transfer stack to generic object
    allocate(rgenericObject%p_cdata(rgenericObject%isize))
    rgenericObject%p_cdata = transfer(rptr, rgenericObject%p_cdata)
    
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(stack_uncast,T)(rgenericObject, p_rstack)

!<description>
    ! This subroutine casts the given generic object into a stack.
!</description>

!<input>
    ! The generic object
    type(t_genericObject), intent(in) :: rgenericObject
!</input>

!<output>
    ! The stack
    type(FEAT2_PP_TEMPLATE_T(t_stack,T)), pointer :: p_rstack
!</output>
!</function>

    ! Internal data structure
    type t_void_ptr
      type(FEAT2_PP_TEMPLATE_TD(t_stack,T,D)), pointer :: p_robj => null()
    end type

    ! Internal variables
    type(t_void_ptr) :: rptr

    if ((rgenericObject%isize .eq. 0) .or.&
        (.not.associated(rgenericObject%p_cdata))) then
      call output_line('Generic object seems to be empty!',&
          OU_CLASS_ERROR,OU_MODE_STD,'stack_uncast')
      call sys_halt()
    end if

    rptr = transfer(rgenericObject%p_cdata, rptr)
    p_rstack => rptr%p_robj

  end subroutine

#endif
