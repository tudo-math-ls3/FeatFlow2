!##############################################################################
!# Tutorial 004h: Datastructures - stacks
!##############################################################################

module tutorial004h

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use stackInt

  implicit none
  private
  
  public :: start_tutorial004h

contains

  ! ***************************************************************************

  subroutine start_tutorial004h

    ! Declare some variables
    type(t_stackInt) :: rstack
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004h")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create stack with space for 10 items
    ! =================================

    call stack_create(rstack, 10)

    ! =================================
    ! Push items on top of stack
    ! =================================

    call stack_push(rstack,  5)
    call stack_push(rstack,  23)
    call stack_push(rstack,  2)
    call stack_push(rstack,  76)
    call stack_push(rstack,  9)
    call stack_push(rstack, -2)
    call stack_push(rstack, -52)
    call stack_push(rstack, -4)

    ! =================================
    ! Check size and content of the stack
    ! =================================
    i = stack_size(rstack)
    call output_line ("Size of the stack: " // trim(sys_siL(i,10)))
    if (stack_empty(rstack)) then
      call output_line ("Stack is empty!")
    else 
      call output_line ("Stack is not empty!")
    end if

    ! =================================
    ! Get the top item from the stack
    ! =================================

    call stack_top(rstack, i)
    call output_line ("Item on top of the stack: " // trim(sys_siL(i,10)))

    ! =================================
    ! Remove (=pop)the top item from the stack
    ! =================================

    call stack_pop(rstack, i)
    call output_line ("Item removed from top of the stack: " // trim(sys_siL(i,10)))

    ! =================================
    ! Get the top item from the stack
    ! =================================

    call stack_top(rstack, i)
    call output_line ("Item on top of the stack: " // trim(sys_siL(i,10)))

    ! =================================
    ! Release stack
    ! =================================

    call stack_release(rstack)

  end subroutine

end module
