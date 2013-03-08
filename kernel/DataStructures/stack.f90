!##############################################################################
!# ****************************************************************************
!# <name> stack </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a meta stack, that is, this data structure provides
!# pointers to all standard stack implementations available.
!#
!# The following routine are available:
!#
!# 1.) stack_init
!#     -> Initialises a meta stack structure
!#
!# 2.) stack_done
!#     -> Finalises a meta stack structure
!#
!# 3.) stack_getbase
!#     -> Returns pointer to a concrete stack
!#
!# </purpose>
!##############################################################################

module stack

!$use omp_lib
  use fsystem
  use genoutput
  use stackDP
  use stackInt
  use stackSP
  use storage

  implicit none

  private
  public :: t_stack
  public :: stack_init
  public :: stack_done
  public :: stack_getbase

  interface stack_getbase
    module procedure stack_getbase_Int
    module procedure stack_getbase_DP
    module procedure stack_getbase_SP
  end interface

!<constants>
!<constantblock description="Global flags for stack implementations">

  ! Stack for integer data
  integer, parameter, public :: STACK_INT    = ST_INT

  ! Stack for double data
  integer, parameter, public :: STACK_DOUBLE = ST_DOUBLE

  ! Stack for single data
  integer, parameter, public :: STACK_SINGLE = ST_SINGLE

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! Meta stack structure
  type t_stack
    private

    ! Pointer to integer-valued stack implementations
    type(t_stackInt),  pointer :: p_StackInt => null()

    ! Pointer to double-valued stack implementations
    type(t_stackDP), pointer :: p_StackDP => null()

    ! Pointer to single-valued stack implementations
    type(t_stackSP), pointer :: p_StackSP => null()

  end type t_stack

!</typeblock>
!</types>

contains

  !************************************************************************

!<subroutine>

  subroutine stack_init(rstack, cstackType)

!<description>
    ! Initialises a meta stack structure
!</description>

!<input>
    ! Stack type
    integer, intent(in) :: cstackType
!</input>

!<output>
    ! Meta stack
    type(t_stack), intent(out) :: rstack
!</output>
!</subroutine>

    select case (cstackType)
    case (STACK_INT)
      allocate(rstack%p_StackInt)

    case (STACK_DOUBLE)
      allocate(rstack%p_StackDP)

    case (STACK_SINGLE)
      allocate(rstack%p_StackSP)

    case DEFAULT
      call output_line('Invalid stack type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'stack_init')
      call sys_halt()
    end select

  end subroutine stack_init

  !************************************************************************

!<subroutine>

  subroutine stack_done(rstack)

!<description>
    ! Finalises a meta stack structure
!</description>

!<inputoutput>
    ! Meta stack
    type(t_stack), intent(inout) :: rstack
!</inputoutput>
!</subroutine>

    if (associated(rstack%p_StackInt)) then
      call stack_release(rstack%p_StackInt)
      deallocate(rstack%p_StackInt)
    end if

    if (associated(rstack%p_StackDP)) then
      call stack_release(rstack%p_StackDP)
      deallocate(rstack%p_StackDP)
    end if

    if (associated(rstack%p_StackSP)) then
      call stack_release(rstack%p_StackSP)
      deallocate(rstack%p_StackSP)
    end if

  end subroutine stack_done

  !************************************************************************

!<subroutine>

  subroutine stack_getbase_Int(rstack, p_rstack)

!<description>
    ! Returns a pointer to the integer-valued stack implementation
!</description>

!<input>
    ! Meta stack
    type(t_stack), intent(in) :: rstack
!</input>

!<output>
    ! Pointer to the stack implementation
    type(t_stackInt), pointer :: p_rstack
!</output>
!</subroutine>

    if (associated(rstack%p_StackInt)) then
      p_rstack => rstack%p_StackInt
    else
      nullify(p_rstack)
    end if

  end subroutine stack_getbase_Int

  !************************************************************************

!<subroutine>

  subroutine stack_getbase_DP(rstack, p_rstack)

!<description>
    ! Returns a pointer to the double-valued stack implementation
!</description>

!<input>
    ! Meta stack
    type(t_stack), intent(in) :: rstack
!</input>

!<output>
    ! Pointer to the stack implementation
    type(t_stackDP), pointer :: p_rstack
!</output>
!</subroutine>

    if (associated(rstack%p_StackDP)) then
      p_rstack => rstack%p_StackDP
    else
      nullify(p_rstack)
    end if

  end subroutine stack_getbase_DP

  !************************************************************************

!<subroutine>

  subroutine stack_getbase_SP(rstack, p_rstack)

!<description>
    ! Returns a pointer to the single-valued stack implementation
!</description>

!<input>
    ! Meta stack
    type(t_stack), intent(in) :: rstack
!</input>

!<output>
    ! Pointer to the stack implementation
    type(t_stackSP), pointer :: p_rstack
!</output>
!</subroutine>

    if (associated(rstack%p_StackSP)) then
      p_rstack => rstack%p_StackSP
    else
      nullify(p_rstack)
    end if

  end subroutine stack_getbase_SP

end module stack
