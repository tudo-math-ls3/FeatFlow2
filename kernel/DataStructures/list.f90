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
  use listInt,       only : t_listInt, list_release
  use listInt_DP,    only : t_listInt_DP, list_release
  use listInt_SP,    only : t_listInt_SP, list_release
  use listInt_Int,   only : t_listInt_Int, list_release
  use listDP,        only : t_listDP, list_release
  use listDP_DP,     only : t_listDP_DP, list_release
  use listDP_SP,     only : t_listDP_SP, list_release
  use listDP_Int,    only : t_listDP_Int, list_release
  use listSP,        only : t_listSP, list_release
  use listSP_DP,     only : t_listSP_DP, list_release
  use listSP_SP,     only : t_listSP_SP, list_release
  use listSP_Int,    only : t_listSP_Int, list_release
  use storage

  implicit none

  private
  public :: t_list
  public :: list_init
  public :: list_done
  public :: list_getbase

  interface list_getbase
    module procedure list_getbase_Int
    module procedure list_getbase_Int_Int
    module procedure list_getbase_Int_DP
    module procedure list_getbase_Int_SP
    module procedure list_getbase_DP
    module procedure list_getbase_DP_Int
    module procedure list_getbase_DP_DP
    module procedure list_getbase_DP_SP
    module procedure list_getbase_SP
    module procedure list_getbase_SP_Int
    module procedure list_getbase_SP_DP
    module procedure list_getbase_SP_SP
  end interface

!<constants>
!<constantblock description="Global flags for list implementations">

  ! list for integer data
  integer, parameter, public :: LIST_INT           = ST_INT
  integer, parameter, public :: LIST_INT_INT       = ST_INT + 100*ST_INT
  integer, parameter, public :: LIST_INT_DOUBLE    = ST_INT + 100*ST_DOUBLE
  integer, parameter, public :: LIST_INT_SINGLE    = ST_INT + 100*ST_SINGLE

  ! list for double data
  integer, parameter, public :: LIST_DOUBLE        = ST_DOUBLE
  integer, parameter, public :: LIST_DOUBLE_INT    = ST_DOUBLE + 100*ST_INT
  integer, parameter, public :: LIST_DOUBLE_DOUBLE = ST_DOUBLE + 100*ST_DOUBLE
  integer, parameter, public :: LIST_DOUBLE_SINGLE = ST_DOUBLE + 100*ST_SINGLE


  ! list for single data
  integer, parameter, public :: LIST_SINGLE        = ST_SINGLE
  integer, parameter, public :: LIST_SINGLE_INT    = ST_SINGLE + 100*ST_INT
  integer, parameter, public :: LIST_SINGLE_DOUBLE = ST_SINGLE + 100*ST_DOUBLE
  integer, parameter, public :: LIST_SINGLE_SINGLE = ST_SINGLE + 100*ST_SINGLE

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! Meta list structure
  type t_list
    private

    ! Pointer to integer-valued list implementations
    type(t_listInt),     pointer :: p_listInt     => null()
    type(t_listInt_Int), pointer :: p_listInt_Int => null()
    type(t_listInt_DP),  pointer :: p_listInt_DP  => null()
    type(t_listInt_SP),  pointer :: p_listInt_SP  => null()

    ! Pointer to double-valued list implementations
    type(t_listDP),      pointer :: p_listDP     => null()
    type(t_listDP_Int),  pointer :: p_listDP_Int => null()
    type(t_listDP_DP),   pointer :: p_listDP_DP  => null()
    type(t_listDP_SP),   pointer :: p_listDP_SP  => null()

    ! Pointer to single-valued list implementations
    type(t_listSP),      pointer :: p_listSP     => null()
    type(t_listSP_Int),  pointer :: p_listSP_Int => null()
    type(t_listSP_DP),   pointer :: p_listSP_DP  => null()
    type(t_listSP_SP),   pointer :: p_listSP_SP  => null()

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
    ! list type
    integer, intent(in) :: clistType
!</input>

!<output>
    ! Meta list
    type(t_list), intent(out) :: rlist
!</output>
!</subroutine>

    select case (clistType)
    case (LIST_INT)
      allocate(rlist%p_listInt)

    case (LIST_INT_INT)
      allocate(rlist%p_listInt_Int)

    case (LIST_INT_DOUBLE)
      allocate(rlist%p_listInt_DP)

    case (LIST_INT_SINGLE)
      allocate(rlist%p_listInt_SP)

    case (LIST_DOUBLE)
      allocate(rlist%p_listDP)

    case (LIST_DOUBLE_INT)
      allocate(rlist%p_listDP_Int)

    case (LIST_DOUBLE_DOUBLE)
      allocate(rlist%p_listDP_DP)

    case (LIST_DOUBLE_SINGLE)
      allocate(rlist%p_listDP_SP)

    case (LIST_SINGLE)
      allocate(rlist%p_listSP)

    case (LIST_SINGLE_INT)
      allocate(rlist%p_listSP_Int)

    case (LIST_SINGLE_DOUBLE)
      allocate(rlist%p_listSP_DP)

    case (LIST_SINGLE_SINGLE)
        allocate(rlist%p_listSP_SP)

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

    if (associated(rlist%p_listInt)) then
      call list_release(rlist%p_listInt)
      deallocate(rlist%p_listInt)
    end if

    if (associated(rlist%p_listInt_Int)) then
      call list_release(rlist%p_listInt_Int)
      deallocate(rlist%p_listInt_INT)
    end if

    if (associated(rlist%p_listInt_DP)) then
      call list_release(rlist%p_listInt_DP)
      deallocate(rlist%p_listInt_DP)
    end if

    if (associated(rlist%p_listInt_SP)) then
      call list_release(rlist%p_listInt_SP)
      deallocate(rlist%p_listInt_SP)
    end if


    if (associated(rlist%p_listDP)) then
      call list_release(rlist%p_listDP)
      deallocate(rlist%p_listDP)
    end if

    if (associated(rlist%p_listDP_Int)) then
      call list_release(rlist%p_listDP_Int)
      deallocate(rlist%p_listDP_Int)
    end if

    if (associated(rlist%p_listDP_DP)) then
      call list_release(rlist%p_listDP_DP)
      deallocate(rlist%p_listDP_DP)
    end if

    if (associated(rlist%p_listDP_SP)) then
      call list_release(rlist%p_listDP_SP)
      deallocate(rlist%p_listDP_SP)
    end if


    if (associated(rlist%p_listSP)) then
      call list_release(rlist%p_listSP)
      deallocate(rlist%p_listSP)
    end if

    if (associated(rlist%p_listSP_Int)) then
      call list_release(rlist%p_listSP_Int)
      deallocate(rlist%p_listSP_Int)
    end if

    if (associated(rlist%p_listSP_DP)) then
      call list_release(rlist%p_listSP_DP)
      deallocate(rlist%p_listSP_DP)
    end if

    if (associated(rlist%p_listSP_SP)) then
      call list_release(rlist%p_listSP_SP)
      deallocate(rlist%p_listSP_SP)
    end if

  end subroutine list_done

  !************************************************************************

!<subroutine>

  subroutine list_getbase_Int(rlist, p_rlist)

!<description>
    ! Returns a pointer to the integer-valued list implementation
    ! without auxiliary data
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

    if (associated(rlist%p_listInt)) then
      p_rlist => rlist%p_listInt
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_Int

  !************************************************************************

!<subroutine>

  subroutine list_getbase_Int_Int(rlist, p_rlist)

!<description>
    ! Returns a pointer to the integer-valued list implementation
    ! with integer-valued auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listInt_Int), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listInt_Int)) then
      p_rlist => rlist%p_listInt_Int
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_Int_Int

  !************************************************************************

!<subroutine>

  subroutine list_getbase_Int_DP(rlist, p_rlist)

!<description>
    ! Returns a pointer to the integer-valued list implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listInt_DP), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listInt_DP)) then
      p_rlist => rlist%p_listInt_DP
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_Int_DP

  !************************************************************************

!<subroutine>

  subroutine list_getbase_Int_SP(rlist, p_rlist)

!<description>
    ! Returns a pointer to the integer-valued list implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listInt_SP), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listInt_SP)) then
      p_rlist => rlist%p_listInt_SP
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_Int_SP

  !************************************************************************

!<subroutine>

  subroutine list_getbase_DP(rlist, p_rlist)

!<description>
    ! Returns a pointer to the double-valued list implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listDP), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listDP)) then
      p_rlist => rlist%p_listDP
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_DP

  !************************************************************************

!<subroutine>

  subroutine list_getbase_DP_Int(rlist, p_rlist)

!<description>
    ! Returns a pointer to the double-valued list implementation
    ! with integer-valued auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listDP_Int), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listDP_Int)) then
      p_rlist => rlist%p_listDP_Int
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_DP_Int

!************************************************************************

!<subroutine>

  subroutine list_getbase_DP_DP(rlist, p_rlist)

!<description>
    ! Returns a pointer to the double-valued list implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listDP_DP), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listDP_DP)) then
      p_rlist => rlist%p_listDP_DP
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_DP_DP

  !************************************************************************

!<subroutine>

  subroutine list_getbase_DP_SP(rlist, p_rlist)

!<description>
    ! Returns a pointer to the double-valued list implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listDP_SP), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listDP_SP)) then
      p_rlist => rlist%p_listDP_SP
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_DP_SP

  !************************************************************************

!<subroutine>

  subroutine list_getbase_SP(rlist, p_rlist)

!<description>
    ! Returns a pointer to the single-valued list implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listSP), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listSP)) then
      p_rlist => rlist%p_listSP
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_SP

  !************************************************************************

!<subroutine>

  subroutine list_getbase_SP_Int(rlist, p_rlist)

!<description>
    ! Returns a pointer to the single-valued list implementation
    ! with integer-valeud auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listSP_Int), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listSP_Int)) then
      p_rlist => rlist%p_listSP_Int
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_SP_Int

  !************************************************************************

!<subroutine>

  subroutine list_getbase_SP_DP(rlist, p_rlist)

!<description>
    ! Returns a pointer to the single-valued list implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listSP_DP), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listSP_DP)) then
      p_rlist => rlist%p_listSP_DP
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_SP_DP

  !************************************************************************

!<subroutine>

  subroutine list_getbase_SP_SP(rlist, p_rlist)

!<description>
    ! Returns a pointer to the single-valued list implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta list
    type(t_list), intent(in) :: rlist
!</input>

!<output>
    ! Pointer to the list implementation
    type(t_listSP_SP), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listSP_SP)) then
      p_rlist => rlist%p_listSP_SP
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_SP_SP

end module list
