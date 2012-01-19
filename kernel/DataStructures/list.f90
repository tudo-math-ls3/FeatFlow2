!##############################################################################
!# ****************************************************************************
!# <name> list </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a meta list, that is, this data structure provides
!# pointers to all standard list implementations available.
!#
!# The following routine are available:
!#
!# 1.) list_init
!#     -> Initialises a meta list structure
!#
!# 2.) list_done
!#     -> Finalises a meta list structure
!#
!# 3.) list_getbase
!#     -> Returns pointer to a concrete list
!#
!# </purpose>
!##############################################################################

module list

!$use omp_lib
  use fsystem
  use genoutput
  use listDble
  use listInt
  use listSngl
  use storage

  implicit none

  interface list_getbase
    module procedure list_getbase_int
    module procedure list_getbase_double
    module procedure list_getbase_single
  end interface

!<constants>
!<constantblock description="Global flags for list implementations">
  
  ! List for integer data
  integer, parameter :: LIST_INT    = ST_INT

  ! List for double data
  integer, parameter :: LIST_DOUBLE = ST_DOUBLE

  ! List for single data
  integer, parameter :: LIST_SINGLE = ST_SINGLE

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! Meta list structure
  type t_list
    private

    ! Pointer to integer-valued list implementations
    type(t_listInt),  pointer :: p_ListInt => null()

    ! Pointer to double-valued list implementations
    type(t_listDble), pointer :: p_ListDble => null()

    ! Pointer to single-valued list implementations
    type(t_listSngl), pointer :: p_ListSngl => null()

  end type t_list

!</typeblock>
!</types>
  
contains

  !************************************************************************

!<subroutine>

  subroutine list_init(rlist, clistType)

!<description>
    ! Initialises a meta list structure
!</description>

!<input>
    ! List type
    integer, intent(in) :: clistType
!</input>

!<output>
    ! Meta list
    type(t_list), intent(out) :: rlist
!</output>
!</subroutine>
    
    select case (clistType)
    case (LIST_INT)
      allocate(rlist%p_ListInt)
      
    case (ST_DOUBLE)
      allocate(rlist%p_ListDble)
    
    case (LIST_SINGLE)
      allocate(rlist%p_ListSngl)
  
    case DEFAULT
      call output_line('Invalid list type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_init')
      call sys_halt()
    end select

  end subroutine list_init

  !************************************************************************

!<subroutine>

  subroutine list_done(rlist)

!<description>
    ! Finalises a meta list structure
!</description>

!<inputoutput>
    ! Meta list
    type(t_list), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    if (associated(rlist%p_ListInt)) then
      call list_release(rlist%p_ListInt)
      deallocate(rlist%p_ListInt)
    end if

    if (associated(rlist%p_ListDble)) then
      call list_release(rlist%p_ListDble)
      deallocate(rlist%p_ListDble)
    end if

    if (associated(rlist%p_ListSngl)) then
      call list_release(rlist%p_ListSngl)
      deallocate(rlist%p_ListSngl)
    end if

  end subroutine list_done

  !************************************************************************

!<subroutine>

  subroutine list_getbase_int(rlist, p_rlist)

!<description>
    ! Returns a pointer to the integer-valued list implementation
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listInt), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_ListInt)) then
      p_rlist => rlist%p_ListInt
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_int

  !************************************************************************

!<subroutine>

  subroutine list_getbase_double(rlist, p_rlist)

!<description>
    ! Returns a pointer to the double-valued list implementation
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listDble), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_ListDble)) then
      p_rlist => rlist%p_ListDble
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_double

  !************************************************************************

!<subroutine>

  subroutine list_getbase_single(rlist, p_rlist)

!<description>
    ! Returns a pointer to the single-valued list implementation
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listSngl), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_ListSngl)) then
      p_rlist => rlist%p_ListSngl
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_single
  
end module list
