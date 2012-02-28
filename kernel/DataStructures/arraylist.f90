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
  use arraylistInt_Dble,  only : t_arraylistInt_Dble, alst_release
  use arraylistInt_Sngl,  only : t_arraylistInt_Sngl, alst_release
  use arraylistInt_Int,   only : t_arraylistInt_Int, alst_release
  use arraylistDble,      only : t_arraylistDble, alst_release
  use arraylistDble_Dble, only : t_arraylistDble_Dble, alst_release
  use arraylistDble_Sngl, only : t_arraylistDble_Sngl, alst_release
  use arraylistDble_Int,  only : t_arraylistDble_Int, alst_release
  use arraylistSngl,      only : t_arraylistSngl, alst_release
  use arraylistSngl_Dble, only : t_arraylistSngl_Dble, alst_release
  use arraylistSngl_Sngl, only : t_arraylistSngl_Sngl, alst_release
  use arraylistSngl_Int,  only : t_arraylistSngl_Int, alst_release
  use storage

  implicit none

  private
  public :: t_arraylist
  public :: alst_init
  public :: alst_done
  public :: alst_getbase

  interface alst_getbase
    module procedure alst_getbase_int
    module procedure alst_getbase_int_int
    module procedure alst_getbase_int_dble
    module procedure alst_getbase_int_sngl
    module procedure alst_getbase_dble
    module procedure alst_getbase_dble_int
    module procedure alst_getbase_dble_dble
    module procedure alst_getbase_dble_sngl
    module procedure alst_getbase_sngl
    module procedure alst_getbase_sngl_int
    module procedure alst_getbase_sngl_dble
    module procedure alst_getbase_sngl_sngl
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
    type(t_arraylistInt_Dble),  pointer :: p_ArraylistInt_Dble  => null()
    type(t_arraylistInt_Sngl),  pointer :: p_ArraylistInt_Sngl  => null()

    ! Pointer to double-valued arraylist implementations
    type(t_arraylistDble),      pointer :: p_ArraylistDble      => null()
    type(t_arraylistDble_Int),  pointer :: p_ArraylistDble_Int  => null()
    type(t_arraylistDble_Dble), pointer :: p_ArraylistDble_Dble => null()
    type(t_arraylistDble_Sngl), pointer :: p_ArraylistDble_Sngl => null()

    ! Pointer to single-valued arraylist implementations
    type(t_arraylistSngl),      pointer :: p_ArraylistSngl      => null()
    type(t_arraylistSngl_Int),  pointer :: p_ArraylistSngl_Int  => null()
    type(t_arraylistSngl_Dble), pointer :: p_ArraylistSngl_Dble => null()
    type(t_arraylistSngl_Sngl), pointer :: p_ArraylistSngl_Sngl => null()

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
      allocate(rarraylist%p_ArraylistInt_Dble)

    case (ALST_INT_SINGLE)
      allocate(rarraylist%p_ArraylistInt_Sngl)
      
    case (ALST_DOUBLE)
      allocate(rarraylist%p_ArraylistDble)

    case (ALST_DOUBLE_INT)
      allocate(rarraylist%p_ArraylistDble_Int)

    case (ALST_DOUBLE_DOUBLE)
      allocate(rarraylist%p_ArraylistDble_Dble)

    case (ALST_DOUBLE_SINGLE)
      allocate(rarraylist%p_ArraylistDble_Sngl)
    
    case (ALST_SINGLE)
      allocate(rarraylist%p_ArraylistSngl)

    case (ALST_SINGLE_INT)
      allocate(rarraylist%p_ArraylistSngl_Int)

    case (ALST_SINGLE_DOUBLE)
      allocate(rarraylist%p_ArraylistSngl_Dble)
      
    case (ALST_SINGLE_SINGLE)
        allocate(rarraylist%p_ArraylistSngl_Sngl)
  
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

    if (associated(rarraylist%p_ArraylistInt_Dble)) then
      call alst_release(rarraylist%p_ArraylistInt_Dble)
      deallocate(rarraylist%p_ArraylistInt_Dble)
    end if

    if (associated(rarraylist%p_ArraylistInt_Sngl)) then
      call alst_release(rarraylist%p_ArraylistInt_Sngl)
      deallocate(rarraylist%p_ArraylistInt_Sngl)
    end if


    if (associated(rarraylist%p_ArraylistDble)) then
      call alst_release(rarraylist%p_ArraylistDble)
      deallocate(rarraylist%p_ArraylistDble)
    end if

    if (associated(rarraylist%p_ArraylistDble_Int)) then
      call alst_release(rarraylist%p_ArraylistDble_Int)
      deallocate(rarraylist%p_ArraylistDble_Int)
    end if

    if (associated(rarraylist%p_ArraylistDble_Dble)) then
      call alst_release(rarraylist%p_ArraylistDble_Dble)
      deallocate(rarraylist%p_ArraylistDble_Dble)
    end if

    if (associated(rarraylist%p_ArraylistDble_Sngl)) then
      call alst_release(rarraylist%p_ArraylistDble_Sngl)
      deallocate(rarraylist%p_ArraylistDble_Sngl)
    end if


    if (associated(rarraylist%p_ArraylistSngl)) then
      call alst_release(rarraylist%p_ArraylistSngl)
      deallocate(rarraylist%p_ArraylistSngl)
    end if

    if (associated(rarraylist%p_ArraylistSngl_Int)) then
      call alst_release(rarraylist%p_ArraylistSngl_Int)
      deallocate(rarraylist%p_ArraylistSngl_Int)
    end if

    if (associated(rarraylist%p_ArraylistSngl_Dble)) then
      call alst_release(rarraylist%p_ArraylistSngl_Dble)
      deallocate(rarraylist%p_ArraylistSngl_Dble)
    end if

    if (associated(rarraylist%p_ArraylistSngl_Sngl)) then
      call alst_release(rarraylist%p_ArraylistSngl_Sngl)
      deallocate(rarraylist%p_ArraylistSngl_Sngl)
    end if

  end subroutine alst_done

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_int(rarraylist, p_rarraylist)

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

  end subroutine alst_getbase_int

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_int_int(rarraylist, p_rarraylist)

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

  end subroutine alst_getbase_int_int

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_int_dble(rarraylist, p_rarraylist)

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
    type(t_arraylistInt_Dble), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistInt_Dble)) then
      p_rarraylist => rarraylist%p_ArraylistInt_Dble
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_int_dble

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_int_sngl(rarraylist, p_rarraylist)

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
    type(t_arraylistInt_Sngl), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistInt_Sngl)) then
      p_rarraylist => rarraylist%p_ArraylistInt_Sngl
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_int_sngl

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_dble(rarraylist, p_rarraylist)

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
    type(t_arraylistDble), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistDble)) then
      p_rarraylist => rarraylist%p_ArraylistDble
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_dble

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_dble_int(rarraylist, p_rarraylist)

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
    type(t_arraylistDble_Int), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistDble_Int)) then
      p_rarraylist => rarraylist%p_ArraylistDble_Int
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_dble_int

!************************************************************************

!<subroutine>

  subroutine alst_getbase_dble_dble(rarraylist, p_rarraylist)

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
    type(t_arraylistDble_Dble), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistDble_Dble)) then
      p_rarraylist => rarraylist%p_ArraylistDble_Dble
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_dble_dble

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_dble_sngl(rarraylist, p_rarraylist)

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
    type(t_arraylistDble_Sngl), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistDble_Sngl)) then
      p_rarraylist => rarraylist%p_ArraylistDble_Sngl
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_dble_sngl

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_sngl(rarraylist, p_rarraylist)

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
    type(t_arraylistSngl), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistSngl)) then
      p_rarraylist => rarraylist%p_ArraylistSngl
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_sngl

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_sngl_int(rarraylist, p_rarraylist)

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
    type(t_arraylistSngl_Int), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistSngl_Int)) then
      p_rarraylist => rarraylist%p_ArraylistSngl_Int
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_sngl_int

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_sngl_dble(rarraylist, p_rarraylist)

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
    type(t_arraylistSngl_Dble), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistSngl_Dble)) then
      p_rarraylist => rarraylist%p_ArraylistSngl_Dble
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_sngl_dble

  !************************************************************************

!<subroutine>

  subroutine alst_getbase_sngl_sngl(rarraylist, p_rarraylist)

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
    type(t_arraylistSngl_Sngl), pointer :: p_rarraylist
!</output>
!</subroutine>

    if (associated(rarraylist%p_ArraylistSngl_Sngl)) then
      p_rarraylist => rarraylist%p_ArraylistSngl_Sngl
    else
      nullify(p_rarraylist)
    end if

  end subroutine alst_getbase_sngl_sngl
  
end module arraylist
