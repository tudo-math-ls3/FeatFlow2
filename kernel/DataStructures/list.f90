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
  use listInt_Dble,  only : t_listInt_Dble, list_release
  use listInt_Sngl,  only : t_listInt_Sngl, list_release
  use listInt_Int,   only : t_listInt_Int, list_release
  use listDble,      only : t_listDble, list_release
  use listDble_Dble, only : t_listDble_Dble, list_release
  use listDble_Sngl, only : t_listDble_Sngl, list_release
  use listDble_Int,  only : t_listDble_Int, list_release
  use listSngl,      only : t_listSngl, list_release
  use listSngl_Dble, only : t_listSngl_Dble, list_release
  use listSngl_Sngl, only : t_listSngl_Sngl, list_release
  use listSngl_Int,  only : t_listSngl_Int, list_release
  use storage

  implicit none

  private
  public :: t_list
  public :: list_init
  public :: list_done
  public :: list_getbase

  interface list_getbase
    module procedure list_getbase_int
    module procedure list_getbase_int_int
    module procedure list_getbase_int_dble
    module procedure list_getbase_int_sngl
    module procedure list_getbase_dble
    module procedure list_getbase_dble_int
    module procedure list_getbase_dble_dble
    module procedure list_getbase_dble_sngl
    module procedure list_getbase_sngl
    module procedure list_getbase_sngl_int
    module procedure list_getbase_sngl_dble
    module procedure list_getbase_sngl_sngl
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
    type(t_listInt),       pointer :: p_listInt       => null()
    type(t_listInt_Int),   pointer :: p_listInt_Int   => null()
    type(t_listInt_Dble),  pointer :: p_listInt_Dble  => null()
    type(t_listInt_Sngl),  pointer :: p_listInt_Sngl  => null()

    ! Pointer to double-valued list implementations
    type(t_listDble),      pointer :: p_listDble      => null()
    type(t_listDble_Int),  pointer :: p_listDble_Int  => null()
    type(t_listDble_Dble), pointer :: p_listDble_Dble => null()
    type(t_listDble_Sngl), pointer :: p_listDble_Sngl => null()

    ! Pointer to single-valued list implementations
    type(t_listSngl),      pointer :: p_listSngl      => null()
    type(t_listSngl_Int),  pointer :: p_listSngl_Int  => null()
    type(t_listSngl_Dble), pointer :: p_listSngl_Dble => null()
    type(t_listSngl_Sngl), pointer :: p_listSngl_Sngl => null()

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
      allocate(rlist%p_listInt_Dble)

    case (LIST_INT_SINGLE)
      allocate(rlist%p_listInt_Sngl)
      
    case (LIST_DOUBLE)
      allocate(rlist%p_listDble)

    case (LIST_DOUBLE_INT)
      allocate(rlist%p_listDble_Int)

    case (LIST_DOUBLE_DOUBLE)
      allocate(rlist%p_listDble_Dble)

    case (LIST_DOUBLE_SINGLE)
      allocate(rlist%p_listDble_Sngl)
    
    case (LIST_SINGLE)
      allocate(rlist%p_listSngl)

    case (LIST_SINGLE_INT)
      allocate(rlist%p_listSngl_Int)

    case (LIST_SINGLE_DOUBLE)
      allocate(rlist%p_listSngl_Dble)
      
    case (LIST_SINGLE_SINGLE)
        allocate(rlist%p_listSngl_Sngl)
  
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

    if (associated(rlist%p_listInt_Dble)) then
      call list_release(rlist%p_listInt_Dble)
      deallocate(rlist%p_listInt_Dble)
    end if

    if (associated(rlist%p_listInt_Sngl)) then
      call list_release(rlist%p_listInt_Sngl)
      deallocate(rlist%p_listInt_Sngl)
    end if


    if (associated(rlist%p_listDble)) then
      call list_release(rlist%p_listDble)
      deallocate(rlist%p_listDble)
    end if

    if (associated(rlist%p_listDble_Int)) then
      call list_release(rlist%p_listDble_Int)
      deallocate(rlist%p_listDble_Int)
    end if

    if (associated(rlist%p_listDble_Dble)) then
      call list_release(rlist%p_listDble_Dble)
      deallocate(rlist%p_listDble_Dble)
    end if

    if (associated(rlist%p_listDble_Sngl)) then
      call list_release(rlist%p_listDble_Sngl)
      deallocate(rlist%p_listDble_Sngl)
    end if


    if (associated(rlist%p_listSngl)) then
      call list_release(rlist%p_listSngl)
      deallocate(rlist%p_listSngl)
    end if

    if (associated(rlist%p_listSngl_Int)) then
      call list_release(rlist%p_listSngl_Int)
      deallocate(rlist%p_listSngl_Int)
    end if

    if (associated(rlist%p_listSngl_Dble)) then
      call list_release(rlist%p_listSngl_Dble)
      deallocate(rlist%p_listSngl_Dble)
    end if

    if (associated(rlist%p_listSngl_Sngl)) then
      call list_release(rlist%p_listSngl_Sngl)
      deallocate(rlist%p_listSngl_Sngl)
    end if

  end subroutine list_done

  !************************************************************************

!<subroutine>

  subroutine list_getbase_int(rlist, p_rlist)

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

  end subroutine list_getbase_int

  !************************************************************************

!<subroutine>

  subroutine list_getbase_int_int(rlist, p_rlist)

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

  end subroutine list_getbase_int_int

  !************************************************************************

!<subroutine>

  subroutine list_getbase_int_dble(rlist, p_rlist)

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
    type(t_listInt_Dble), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listInt_Dble)) then
      p_rlist => rlist%p_listInt_Dble
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_int_dble

  !************************************************************************

!<subroutine>

  subroutine list_getbase_int_sngl(rlist, p_rlist)

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
    type(t_listInt_Sngl), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listInt_Sngl)) then
      p_rlist => rlist%p_listInt_Sngl
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_int_sngl

  !************************************************************************

!<subroutine>

  subroutine list_getbase_dble(rlist, p_rlist)

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
    type(t_listDble), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listDble)) then
      p_rlist => rlist%p_listDble
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_dble

  !************************************************************************

!<subroutine>

  subroutine list_getbase_dble_int(rlist, p_rlist)

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
    type(t_listDble_Int), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listDble_Int)) then
      p_rlist => rlist%p_listDble_Int
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_dble_int

!************************************************************************

!<subroutine>

  subroutine list_getbase_dble_dble(rlist, p_rlist)

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
    type(t_listDble_Dble), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listDble_Dble)) then
      p_rlist => rlist%p_listDble_Dble
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_dble_dble

  !************************************************************************

!<subroutine>

  subroutine list_getbase_dble_sngl(rlist, p_rlist)

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
    type(t_listDble_Sngl), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listDble_Sngl)) then
      p_rlist => rlist%p_listDble_Sngl
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_dble_sngl

  !************************************************************************

!<subroutine>

  subroutine list_getbase_sngl(rlist, p_rlist)

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
    type(t_listSngl), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listSngl)) then
      p_rlist => rlist%p_listSngl
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_sngl

  !************************************************************************

!<subroutine>

  subroutine list_getbase_sngl_int(rlist, p_rlist)

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
    type(t_listSngl_Int), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listSngl_Int)) then
      p_rlist => rlist%p_listSngl_Int
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_sngl_int

  !************************************************************************

!<subroutine>

  subroutine list_getbase_sngl_dble(rlist, p_rlist)

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
    type(t_listSngl_Dble), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listSngl_Dble)) then
      p_rlist => rlist%p_listSngl_Dble
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_sngl_dble

  !************************************************************************

!<subroutine>

  subroutine list_getbase_sngl_sngl(rlist, p_rlist)

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
    type(t_listSngl_Sngl), pointer :: p_rlist
!</output>
!</subroutine>

    if (associated(rlist%p_listSngl_Sngl)) then
      p_rlist => rlist%p_listSngl_Sngl
    else
      nullify(p_rlist)
    end if

  end subroutine list_getbase_sngl_sngl
  
end module list
