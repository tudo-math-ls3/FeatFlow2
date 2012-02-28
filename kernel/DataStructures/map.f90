!##############################################################################
!# ****************************************************************************
!# <name> map </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a meta map, that is, this data structure provides
!# pointers to all standard map implementations available.
!#
!# The following routine are available:
!#
!# 1.) map_init
!#     -> Initialises a meta map structure
!#
!# 2.) map_done
!#     -> Finalises a meta map structure
!#
!# 3.) map_getbase
!#     -> Returns pointer to a concrete map
!#
!# </purpose>
!##############################################################################

module map

!$use omp_lib
  use fsystem
  use genoutput
  use mapInt
  use mapInt_Dble,  only : t_mapInt_Dble, map_release
  use mapInt_Sngl,  only : t_mapInt_Sngl, map_release
  use mapInt_Int,   only : t_mapInt_Int, map_release
  use mapDble,      only : t_mapDble, map_release
  use mapDble_Dble, only : t_mapDble_Dble, map_release
  use mapDble_Sngl, only : t_mapDble_Sngl, map_release
  use mapDble_Int,  only : t_mapDble_Int, map_release
  use mapSngl,      only : t_mapSngl, map_release
  use mapSngl_Dble, only : t_mapSngl_Dble, map_release
  use mapSngl_Sngl, only : t_mapSngl_Sngl, map_release
  use mapSngl_Int,  only : t_mapSngl_Int, map_release
  use storage

  implicit none

  private
  public :: t_map
  public :: map_init
  public :: map_done
  public :: map_getbase

  interface map_getbase
    module procedure map_getbase_int
    module procedure map_getbase_int_int
    module procedure map_getbase_int_dble
    module procedure map_getbase_int_sngl
    module procedure map_getbase_dble
    module procedure map_getbase_dble_int
    module procedure map_getbase_dble_dble
    module procedure map_getbase_dble_sngl
    module procedure map_getbase_sngl
    module procedure map_getbase_sngl_int
    module procedure map_getbase_sngl_dble
    module procedure map_getbase_sngl_sngl
  end interface

!<constants>
!<constantblock description="Global flags for map implementations">
  
  ! map for integer data
  integer, parameter, public :: MAP_INT           = ST_INT
  integer, parameter, public :: MAP_INT_INT       = ST_INT + 100*ST_INT
  integer, parameter, public :: MAP_INT_DOUBLE    = ST_INT + 100*ST_DOUBLE
  integer, parameter, public :: MAP_INT_SINGLE    = ST_INT + 100*ST_SINGLE

  ! map for double data
  integer, parameter, public :: MAP_DOUBLE        = ST_DOUBLE
  integer, parameter, public :: MAP_DOUBLE_INT    = ST_DOUBLE + 100*ST_INT
  integer, parameter, public :: MAP_DOUBLE_DOUBLE = ST_DOUBLE + 100*ST_DOUBLE
  integer, parameter, public :: MAP_DOUBLE_SINGLE = ST_DOUBLE + 100*ST_SINGLE


  ! map for single data
  integer, parameter, public :: MAP_SINGLE        = ST_SINGLE
  integer, parameter, public :: MAP_SINGLE_INT    = ST_SINGLE + 100*ST_INT
  integer, parameter, public :: MAP_SINGLE_DOUBLE = ST_SINGLE + 100*ST_DOUBLE
  integer, parameter, public :: MAP_SINGLE_SINGLE = ST_SINGLE + 100*ST_SINGLE

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! Meta map structure
  type t_map
    private

    ! Pointer to integer-valued map implementations
    type(t_mapInt),       pointer :: p_mapInt       => null()
    type(t_mapInt_Int),   pointer :: p_mapInt_Int   => null()
    type(t_mapInt_Dble),  pointer :: p_mapInt_Dble  => null()
    type(t_mapInt_Sngl),  pointer :: p_mapInt_Sngl  => null()

    ! Pointer to double-valued map implementations
    type(t_mapDble),      pointer :: p_mapDble      => null()
    type(t_mapDble_Int),  pointer :: p_mapDble_Int  => null()
    type(t_mapDble_Dble), pointer :: p_mapDble_Dble => null()
    type(t_mapDble_Sngl), pointer :: p_mapDble_Sngl => null()

    ! Pointer to single-valued map implementations
    type(t_mapSngl),      pointer :: p_mapSngl      => null()
    type(t_mapSngl_Int),  pointer :: p_mapSngl_Int  => null()
    type(t_mapSngl_Dble), pointer :: p_mapSngl_Dble => null()
    type(t_mapSngl_Sngl), pointer :: p_mapSngl_Sngl => null()

  end type t_map

!</typeblock>
!</types>
  
contains

  !************************************************************************

!<subroutine>

  subroutine map_init(rmap, cmapType)

!<description>
    ! Initialises a meta map structure
!</description>

!<input>
    ! map type
    integer, intent(in) :: cmapType
!</input>

!<output>
    ! Meta map
    type(t_map), intent(out) :: rmap
!</output>
!</subroutine>
    
    select case (cmapType)
    case (MAP_INT)
      allocate(rmap%p_mapInt)

    case (MAP_INT_INT)
      allocate(rmap%p_mapInt_Int)

    case (MAP_INT_DOUBLE)
      allocate(rmap%p_mapInt_Dble)

    case (MAP_INT_SINGLE)
      allocate(rmap%p_mapInt_Sngl)
      
    case (MAP_DOUBLE)
      allocate(rmap%p_mapDble)

    case (MAP_DOUBLE_INT)
      allocate(rmap%p_mapDble_Int)

    case (MAP_DOUBLE_DOUBLE)
      allocate(rmap%p_mapDble_Dble)

    case (MAP_DOUBLE_SINGLE)
      allocate(rmap%p_mapDble_Sngl)
    
    case (MAP_SINGLE)
      allocate(rmap%p_mapSngl)

    case (MAP_SINGLE_INT)
      allocate(rmap%p_mapSngl_Int)

    case (MAP_SINGLE_DOUBLE)
      allocate(rmap%p_mapSngl_Dble)
      
    case (MAP_SINGLE_SINGLE)
        allocate(rmap%p_mapSngl_Sngl)
  
    case DEFAULT
      call output_line('Invalid map type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'map_init')
      call sys_halt()
    end select

  end subroutine map_init

  !************************************************************************

!<subroutine>

  subroutine map_done(rmap)

!<description>
    ! Finalises a meta map structure
!</description>

!<inputoutput>
    ! Meta map
    type(t_map), intent(inout) :: rmap
!</inputoutput>
!</subroutine>

    if (associated(rmap%p_mapInt)) then
      call map_release(rmap%p_mapInt)
      deallocate(rmap%p_mapInt)
    end if

    if (associated(rmap%p_mapInt_Int)) then
      call map_release(rmap%p_mapInt_Int)
      deallocate(rmap%p_mapInt_INT)
    end if

    if (associated(rmap%p_mapInt_Dble)) then
      call map_release(rmap%p_mapInt_Dble)
      deallocate(rmap%p_mapInt_Dble)
    end if

    if (associated(rmap%p_mapInt_Sngl)) then
      call map_release(rmap%p_mapInt_Sngl)
      deallocate(rmap%p_mapInt_Sngl)
    end if


    if (associated(rmap%p_mapDble)) then
      call map_release(rmap%p_mapDble)
      deallocate(rmap%p_mapDble)
    end if

    if (associated(rmap%p_mapDble_Int)) then
      call map_release(rmap%p_mapDble_Int)
      deallocate(rmap%p_mapDble_Int)
    end if

    if (associated(rmap%p_mapDble_Dble)) then
      call map_release(rmap%p_mapDble_Dble)
      deallocate(rmap%p_mapDble_Dble)
    end if

    if (associated(rmap%p_mapDble_Sngl)) then
      call map_release(rmap%p_mapDble_Sngl)
      deallocate(rmap%p_mapDble_Sngl)
    end if


    if (associated(rmap%p_mapSngl)) then
      call map_release(rmap%p_mapSngl)
      deallocate(rmap%p_mapSngl)
    end if

    if (associated(rmap%p_mapSngl_Int)) then
      call map_release(rmap%p_mapSngl_Int)
      deallocate(rmap%p_mapSngl_Int)
    end if

    if (associated(rmap%p_mapSngl_Dble)) then
      call map_release(rmap%p_mapSngl_Dble)
      deallocate(rmap%p_mapSngl_Dble)
    end if

    if (associated(rmap%p_mapSngl_Sngl)) then
      call map_release(rmap%p_mapSngl_Sngl)
      deallocate(rmap%p_mapSngl_Sngl)
    end if

  end subroutine map_done

  !************************************************************************

!<subroutine>

  subroutine map_getbase_int(rmap, p_rmap)

!<description>
    ! Returns a pointer to the integer-valued map implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapInt), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapInt)) then
      p_rmap => rmap%p_mapInt
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_int

  !************************************************************************

!<subroutine>

  subroutine map_getbase_int_int(rmap, p_rmap)

!<description>
    ! Returns a pointer to the integer-valued map implementation
    ! with integer-valued auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapInt_Int), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapInt_Int)) then
      p_rmap => rmap%p_mapInt_Int
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_int_int

  !************************************************************************

!<subroutine>

  subroutine map_getbase_int_dble(rmap, p_rmap)

!<description>
    ! Returns a pointer to the integer-valued map implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapInt_Dble), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapInt_Dble)) then
      p_rmap => rmap%p_mapInt_Dble
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_int_dble

  !************************************************************************

!<subroutine>

  subroutine map_getbase_int_sngl(rmap, p_rmap)

!<description>
    ! Returns a pointer to the integer-valued map implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapInt_Sngl), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapInt_Sngl)) then
      p_rmap => rmap%p_mapInt_Sngl
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_int_sngl

  !************************************************************************

!<subroutine>

  subroutine map_getbase_dble(rmap, p_rmap)

!<description>
    ! Returns a pointer to the double-valued map implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapDble), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapDble)) then
      p_rmap => rmap%p_mapDble
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_dble

  !************************************************************************

!<subroutine>

  subroutine map_getbase_dble_int(rmap, p_rmap)

!<description>
    ! Returns a pointer to the double-valued map implementation
    ! with integer-valued auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapDble_Int), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapDble_Int)) then
      p_rmap => rmap%p_mapDble_Int
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_dble_int

!************************************************************************

!<subroutine>

  subroutine map_getbase_dble_dble(rmap, p_rmap)

!<description>
    ! Returns a pointer to the double-valued map implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapDble_Dble), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapDble_Dble)) then
      p_rmap => rmap%p_mapDble_Dble
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_dble_dble

  !************************************************************************

!<subroutine>

  subroutine map_getbase_dble_sngl(rmap, p_rmap)

!<description>
    ! Returns a pointer to the double-valued map implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapDble_Sngl), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapDble_Sngl)) then
      p_rmap => rmap%p_mapDble_Sngl
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_dble_sngl

  !************************************************************************

!<subroutine>

  subroutine map_getbase_sngl(rmap, p_rmap)

!<description>
    ! Returns a pointer to the single-valued map implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapSngl), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapSngl)) then
      p_rmap => rmap%p_mapSngl
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_sngl

  !************************************************************************

!<subroutine>

  subroutine map_getbase_sngl_int(rmap, p_rmap)

!<description>
    ! Returns a pointer to the single-valued map implementation
    ! with integer-valeud auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapSngl_Int), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapSngl_Int)) then
      p_rmap => rmap%p_mapSngl_Int
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_sngl_int

  !************************************************************************

!<subroutine>

  subroutine map_getbase_sngl_dble(rmap, p_rmap)

!<description>
    ! Returns a pointer to the single-valued map implementation
    ! with double-valued auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapSngl_Dble), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapSngl_Dble)) then
      p_rmap => rmap%p_mapSngl_Dble
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_sngl_dble

  !************************************************************************

!<subroutine>

  subroutine map_getbase_sngl_sngl(rmap, p_rmap)

!<description>
    ! Returns a pointer to the single-valued map implementation
    ! with single-valued auxiliary data
!</description>

!<input>
    ! Meta map
    type(t_map), intent(in) :: rmap
!</input>

!<output>
    ! Pointer to the map implementation
    type(t_mapSngl_Sngl), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapSngl_Sngl)) then
      p_rmap => rmap%p_mapSngl_Sngl
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_sngl_sngl
  
end module map
