!##############################################################################
!# ****************************************************************************
!# <name> arraylist </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a meta arraylist, that is, this data structure provides
!# pointers to all standard arraylist implementations available.
!#
!# The following routine are available:
!#
!# 1.) alst_init
!#     -> Initialises a meta arraylist structure
!#
!# 2.) alst_done
!#     -> Finalises a meta arraylist structure
!#
!# 3.) alst_getbase
!#     -> Returns pointer to a concrete arraylist
!#
!# </purpose>
!##############################################################################

module arraylist

!$use omp_lib
  use fsystem
  use genoutput
  use arraylistInt,       only : t_arraylistInt, alst_release
  use arraylistInt_DP,    only : t_arraylistInt_DP, alst_release
  use arraylistInt_SP,    only : t_arraylistInt_SP, alst_release
  use arraylistInt_Int,   only : t_arraylistInt_Int, alst_release
  use arraylistDP,        only : t_arraylistDP, alst_release
  use arraylistDP_DP,     only : t_arraylistDP_DP, alst_release
  use arraylistDP_SP,     only : t_arraylistDP_SP, alst_release
  use arraylistDP_Int,    only : t_arraylistDP_Int, alst_release
  use arraylistSP,        only : t_arraylistSP, alst_release
  use arraylistSP_DP,     only : t_arraylistSP_DP, alst_release
  use arraylistSP_SP,     only : t_arraylistSP_SP, alst_release
  use arraylistSP_Int,    only : t_arraylistSP_Int, alst_release
  use storage

  implicit none

  private
  public :: t_arraylist
  public :: alst_init
  public :: alst_done
  public :: alst_getbase

  interface alst_getbase
    module procedure alst_getbase_Int
    module procedure alst_getbase_Int_Int
    module procedure alst_getbase_Int_DP
    module procedure alst_getbase_Int_SP
    module procedure alst_getbase_DP
    module procedure alst_getbase_DP_Int
    module procedure alst_getbase_DP_DP
    module procedure alst_getbase_DP_SP
    module procedure alst_getbase_SP
    module procedure alst_getbase_SP_Int
    module procedure alst_getbase_SP_DP
    module procedure alst_getbase_SP_SP
  end interface

!<constants>
!<constantblock description="Global flags for arraylist implementations">

  ! Arraylist for integer data
  integer, parameter :: ALST_INT           = ST_INT
  integer, parameter :: ALST_INT_INT       = ST_INT + 100*ST_INT
  integer, parameter :: ALST_INT_DOUBLE    = ST_INT + 100*ST_DOUBLE
  integer, parameter :: ALST_INT_SINGLE    = ST_INT + 100*ST_SINGLE

  ! Arraylist for double data
  integer, parameter :: ALST_DOUBLE        = ST_DOUBLE
  integer, parameter :: ALST_DOUBLE_INT    = ST_DOUBLE + 100*ST_INT
  integer, parameter :: ALST_DOUBLE_DOUBLE = ST_DOUBLE + 100*ST_DOUBLE
  integer, parameter :: ALST_DOUBLE_SINGLE = ST_DOUBLE + 100*ST_SINGLE


  ! Arraylist for single data
  integer, parameter :: ALST_SINGLE        = ST_SINGLE
  integer, parameter :: ALST_SINGLE_INT    = ST_SINGLE + 100*ST_INT
  integer, parameter :: ALST_SINGLE_DOUBLE = ST_SINGLE + 100*ST_DOUBLE
  integer, parameter :: ALST_SINGLE_SINGLE = ST_SINGLE + 100*ST_SINGLE

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! Meta arraylist structure
  type t_arraylist
    private

    ! Pointer to integer-valued arraylist implementations
    type(t_arraylistInt),       pointer :: p_ArraylistInt       => null()
    type(t_arraylistInt_Int),   pointer :: p_ArraylistInt_Int   => null()
    type(t_arraylistInt_DP),  pointer :: p_ArraylistInt_DP  => null()
    type(t_arraylistInt_SP),  pointer :: p_ArraylistInt_SP  => null()

    ! Pointer to double-valued arraylist implementations
    type(t_arraylistDP),      pointer :: p_ArraylistDP      => null()
    type(t_arraylistDP_Int),  pointer :: p_ArraylistDP_Int  => null()
    type(t_arraylistDP_DP), pointer :: p_ArraylistDP_DP => null()
    type(t_arraylistDP_SP), pointer :: p_ArraylistDP_SP => null()

    ! Pointer to single-valued arraylist implementations
    type(t_arraylistSP),      pointer :: p_ArraylistSP      => null()
    type(t_arraylistSP_Int),  pointer :: p_ArraylistSP_Int  => null()
    type(t_arraylistSP_DP), pointer :: p_ArraylistSP_DP => null()
    type(t_arraylistSP_SP), pointer :: p_ArraylistSP_SP => null()

  end type t_arraylist

!</typeblock>
!</types>

contains

  !************************************************************************

!<subroutine>

  subroutine alst_init(rarraylist, carraylistType)

!<description>
    ! Initialises a meta arraylist structure
!</description>

!<input>
    ! Arraylist type
    integer, intent(in) :: carraylistType
!</input>

!<output>
    ! Meta arraylist
    type(t_arraylist), intent(out) :: rarraylist
!</output>
!</subroutine>

    select case (carraylistType)
    case (ALST_INT)
      allocate(rarraylist%p_ArraylistInt)

    case (ALST_INT_INT)
      allocate(rarraylist%p_ArraylistInt_Int)

    case (ALST_INT_DOUBLE)
      allocate(rarraylist%p_ArraylistInt_DP)

    case (ALST_INT_SINGLE)
      allocate(rarraylist%p_ArraylistInt_SP)

    case (ALST_DOUBLE)
      allocate(rarraylist%p_ArraylistDP)

    case (ALST_DOUBLE_INT)
      allocate(rarraylist%p_ArraylistDP_Int)

    case (ALST_DOUBLE_DOUBLE)
      allocate(rarraylist%p_ArraylistDP_DP)

    case (ALST_DOUBLE_SINGLE)
      allocate(rarraylist%p_ArraylistDP_SP)

    case (ALST_SINGLE)
      allocate(rarraylist%p_ArraylistSP)

    case (ALST_SINGLE_INT)
      allocate(rarraylist%p_ArraylistSP_Int)

    case (ALST_SINGLE_DOUBLE)
      allocate(rarraylist%p_ArraylistSP_DP)

    case (ALST_SINGLE_SINGLE)
        allocate(rarraylist%p_ArraylistSP_SP)

    case DEFAULT
      call output_line('Invalid arraylist type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'alst_init')
      call sys_halt()
    end select

  end subroutine alst_init

  !************************************************************************

!<subroutine>

  subroutine alst_done(rarraylist)

!<description>
    ! Finalises a meta arraylist structure
!</description>

!<inputoutput>
    ! Meta arraylist
    type(t_arraylist), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    if (associated(rarraylist%p_ArraylistInt)) then
      call alst_release(rarraylist%p_ArraylistInt)
      deallocate(rarraylist%p_ArraylistInt)
    end if

    if (associated(rarraylist%p_ArraylistInt_Int)) then
      call alst_release(rarraylist%p_ArraylistInt_Int)
      deallocate(rarraylist%p_ArraylistInt_INT)
    end if

    if (associated(rarraylist%p_ArraylistInt_DP)) then
      call alst_release(rarraylist%p_ArraylistInt_DP)
      deallocate(rarraylist%p_ArraylistInt_DP)
    end if

    if (associated(rarraylist%p_ArraylistInt_SP)) then
      call alst_release(rarraylist%p_ArraylistInt_SP)
      deallocate(rarraylist%p_ArraylistInt_SP)
    end if


    if (associated(rarraylist%p_ArraylistDP)) then
      call alst_release(rarraylist%p_ArraylistDP)
      deallocate(rarraylist%p_ArraylistDP)
    end if

    if (associated(rarraylist%p_ArraylistDP_Int)) then
      call alst_release(rarraylist%p_ArraylistDP_Int)
      deallocate(rarraylist%p_ArraylistDP_Int)
    end if

    if (associated(rarraylist%p_ArraylistDP_DP)) then
      call alst_release(rarraylist%p_ArraylistDP_DP)
      deallocate(rarraylist%p_ArraylistDP_DP)
    end if

    if (associated(rarraylist%p_ArraylistDP_SP)) then
      call alst_release(rarraylist%p_ArraylistDP_SP)
      deallocate(rarraylist%p_ArraylistDP_SP)
    end if


    if (associated(rarraylist%p_ArraylistSP)) then
      call alst_release(rarraylist%p_ArraylistSP)
      deallocate(rarraylist%p_ArraylistSP)
    end if

    if (associated(rarraylist%p_ArraylistSP_Int)) then
      call alst_release(rarraylist%p_ArraylistSP_Int)
      deallocate(rarraylist%p_ArraylistSP_Int)
    end if

    if (associated(rarraylist%p_ArraylistSP_DP)) then
      call alst_release(rarraylist%p_ArraylistSP_DP)
      deallocate(rarraylist%p_ArraylistSP_DP)
    end if

    if (associated(rarraylist%p_ArraylistSP_SP)) then
      call alst_release(rarraylist%p_ArraylistSP_SP)
      deallocate(rarraylist%p_ArraylistSP_SP)
    end if

  end subroutine alst_done

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_Int(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the integer-valued arraylist implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistInt), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistInt)) then
      p_rarraylist => rarraylist%p_ArraylistInt
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_Int

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_Int_Int(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the integer-valued arraylist implementation
    ! with integer-valued auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistInt_Int), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistInt_Int)) then
      p_rarraylist => rarraylist%p_ArraylistInt_Int
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_Int_Int

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_Int_DP(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the integer-valued arraylist implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistInt_DP), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistInt_DP)) then
      p_rarraylist => rarraylist%p_ArraylistInt_DP
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_Int_DP

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_Int_SP(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the integer-valued arraylist implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistInt_SP), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistInt_SP)) then
      p_rarraylist => rarraylist%p_ArraylistInt_SP
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_Int_SP

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_DP(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the double-valued arraylist implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistDP), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistDP)) then
      p_rarraylist => rarraylist%p_ArraylistDP
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_DP

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_DP_Int(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the double-valued arraylist implementation
    ! with integer-valued auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistDP_Int), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistDP_Int)) then
      p_rarraylist => rarraylist%p_ArraylistDP_Int
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_DP_Int

!************************************************************************

!<subroutine>

  subroutine alst_getbase_DP_DP(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the double-valued arraylist implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistDP_DP), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistDP_DP)) then
      p_rarraylist => rarraylist%p_ArraylistDP_DP
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_DP_DP

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_DP_SP(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the double-valued arraylist implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistDP_SP), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistDP_SP)) then
      p_rarraylist => rarraylist%p_ArraylistDP_SP
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_DP_SP

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_SP(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the single-valued arraylist implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistSP), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistSP)) then
      p_rarraylist => rarraylist%p_ArraylistSP
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_SP

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_SP_Int(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the single-valued arraylist implementation
    ! with integer-valeud auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistSP_Int), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistSP_Int)) then
      p_rarraylist => rarraylist%p_ArraylistSP_Int
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_SP_Int

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_SP_DP(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the single-valued arraylist implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistSP_DP), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistSP_DP)) then
      p_rarraylist => rarraylist%p_ArraylistSP_DP
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_SP_DP

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_SP_SP(rarraylist, p_rarraylist)

!<description>
    ! Returns a pointer to the single-valued arraylist implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta arraylist
    type(t_arraylist), intent(in) :: rarraylist
!</input>

!<output>
    ! Pointer to the arraylist implementation
    type(t_arraylistSP_SP), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistSP_SP)) then
      p_rarraylist => rarraylist%p_ArraylistSP_SP
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_SP_SP

end module arraylist
