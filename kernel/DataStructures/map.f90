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
  use mapInt_DP,    only : t_mapInt_DP, map_release
  use mapInt_SP,    only : t_mapInt_SP, map_release
  use mapInt_Int,   only : t_mapInt_Int, map_release
  use mapDP,        only : t_mapDP, map_release
  use mapDP_DP,     only : t_mapDP_DP, map_release
  use mapDP_SP,     only : t_mapDP_SP, map_release
  use mapDP_Int,    only : t_mapDP_Int, map_release
  use mapSP,        only : t_mapSP, map_release
  use mapSP_DP,     only : t_mapSP_DP, map_release
  use mapSP_SP,     only : t_mapSP_SP, map_release
  use mapSP_Int,    only : t_mapSP_Int, map_release
  use storage

  implicit none

  private
  public :: t_map
  public :: map_init
  public :: map_done
  public :: map_getbase

  interface map_getbase
    module procedure map_getbase_Int
    module procedure map_getbase_Int_Int
    module procedure map_getbase_Int_DP
    module procedure map_getbase_Int_SP
    module procedure map_getbase_DP
    module procedure map_getbase_DP_Int
    module procedure map_getbase_DP_DP
    module procedure map_getbase_DP_SP
    module procedure map_getbase_SP
    module procedure map_getbase_SP_Int
    module procedure map_getbase_SP_DP
    module procedure map_getbase_SP_SP
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
    type(t_mapInt_DP),  pointer :: p_mapInt_DP  => null()
    type(t_mapInt_SP),  pointer :: p_mapInt_SP  => null()

    ! Pointer to double-valued map implementations
    type(t_mapDP),      pointer :: p_mapDP      => null()
    type(t_mapDP_Int),  pointer :: p_mapDP_Int  => null()
    type(t_mapDP_DP), pointer :: p_mapDP_DP => null()
    type(t_mapDP_SP), pointer :: p_mapDP_SP => null()

    ! Pointer to single-valued map implementations
    type(t_mapSP),      pointer :: p_mapSP      => null()
    type(t_mapSP_Int),  pointer :: p_mapSP_Int  => null()
    type(t_mapSP_DP), pointer :: p_mapSP_DP => null()
    type(t_mapSP_SP), pointer :: p_mapSP_SP => null()

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
      allocate(rmap%p_mapInt_DP)

    case (MAP_INT_SINGLE)
      allocate(rmap%p_mapInt_SP)

    case (MAP_DOUBLE)
      allocate(rmap%p_mapDP)

    case (MAP_DOUBLE_INT)
      allocate(rmap%p_mapDP_Int)

    case (MAP_DOUBLE_DOUBLE)
      allocate(rmap%p_mapDP_DP)

    case (MAP_DOUBLE_SINGLE)
      allocate(rmap%p_mapDP_SP)

    case (MAP_SINGLE)
      allocate(rmap%p_mapSP)

    case (MAP_SINGLE_INT)
      allocate(rmap%p_mapSP_Int)

    case (MAP_SINGLE_DOUBLE)
      allocate(rmap%p_mapSP_DP)

    case (MAP_SINGLE_SINGLE)
        allocate(rmap%p_mapSP_SP)

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

    if (associated(rmap%p_mapInt_DP)) then
      call map_release(rmap%p_mapInt_DP)
      deallocate(rmap%p_mapInt_DP)
    end if

    if (associated(rmap%p_mapInt_SP)) then
      call map_release(rmap%p_mapInt_SP)
      deallocate(rmap%p_mapInt_SP)
    end if


    if (associated(rmap%p_mapDP)) then
      call map_release(rmap%p_mapDP)
      deallocate(rmap%p_mapDP)
    end if

    if (associated(rmap%p_mapDP_Int)) then
      call map_release(rmap%p_mapDP_Int)
      deallocate(rmap%p_mapDP_Int)
    end if

    if (associated(rmap%p_mapDP_DP)) then
      call map_release(rmap%p_mapDP_DP)
      deallocate(rmap%p_mapDP_DP)
    end if

    if (associated(rmap%p_mapDP_SP)) then
      call map_release(rmap%p_mapDP_SP)
      deallocate(rmap%p_mapDP_SP)
    end if


    if (associated(rmap%p_mapSP)) then
      call map_release(rmap%p_mapSP)
      deallocate(rmap%p_mapSP)
    end if

    if (associated(rmap%p_mapSP_Int)) then
      call map_release(rmap%p_mapSP_Int)
      deallocate(rmap%p_mapSP_Int)
    end if

    if (associated(rmap%p_mapSP_DP)) then
      call map_release(rmap%p_mapSP_DP)
      deallocate(rmap%p_mapSP_DP)
    end if

    if (associated(rmap%p_mapSP_SP)) then
      call map_release(rmap%p_mapSP_SP)
      deallocate(rmap%p_mapSP_SP)
    end if

  end subroutine map_done

  !************************************************************************

!<subroutine>

  subroutine map_getbase_Int(rmap, p_rmap)

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

  end subroutine map_getbase_Int

  !************************************************************************

!<subroutine>

  subroutine map_getbase_Int_Int(rmap, p_rmap)

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

  end subroutine map_getbase_Int_Int

  !************************************************************************

!<subroutine>

  subroutine map_getbase_Int_DP(rmap, p_rmap)

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
    type(t_mapInt_DP), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapInt_DP)) then
      p_rmap => rmap%p_mapInt_DP
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_Int_DP

  !************************************************************************

!<subroutine>

  subroutine map_getbase_Int_SP(rmap, p_rmap)

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
    type(t_mapInt_SP), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapInt_SP)) then
      p_rmap => rmap%p_mapInt_SP
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_Int_SP

  !************************************************************************

!<subroutine>

  subroutine map_getbase_DP(rmap, p_rmap)

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
    type(t_mapDP), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapDP)) then
      p_rmap => rmap%p_mapDP
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_DP

  !************************************************************************

!<subroutine>

  subroutine map_getbase_DP_Int(rmap, p_rmap)

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
    type(t_mapDP_Int), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapDP_Int)) then
      p_rmap => rmap%p_mapDP_Int
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_DP_Int

!************************************************************************

!<subroutine>

  subroutine map_getbase_DP_DP(rmap, p_rmap)

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
    type(t_mapDP_DP), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapDP_DP)) then
      p_rmap => rmap%p_mapDP_DP
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_DP_DP

  !************************************************************************

!<subroutine>

  subroutine map_getbase_DP_SP(rmap, p_rmap)

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
    type(t_mapDP_SP), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapDP_SP)) then
      p_rmap => rmap%p_mapDP_SP
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_DP_SP

  !************************************************************************

!<subroutine>

  subroutine map_getbase_SP(rmap, p_rmap)

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
    type(t_mapSP), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapSP)) then
      p_rmap => rmap%p_mapSP
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_SP

  !************************************************************************

!<subroutine>

  subroutine map_getbase_SP_Int(rmap, p_rmap)

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
    type(t_mapSP_Int), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapSP_Int)) then
      p_rmap => rmap%p_mapSP_Int
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_SP_Int

  !************************************************************************

!<subroutine>

  subroutine map_getbase_SP_DP(rmap, p_rmap)

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
    type(t_mapSP_DP), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapSP_DP)) then
      p_rmap => rmap%p_mapSP_DP
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_SP_DP

  !************************************************************************

!<subroutine>

  subroutine map_getbase_SP_SP(rmap, p_rmap)

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
    type(t_mapSP_SP), pointer :: p_rmap
!</output>
!</subroutine>

    if (associated(rmap%p_mapSP_SP)) then
      p_rmap => rmap%p_mapSP_SP
    else
      nullify(p_rmap)
    end if

  end subroutine map_getbase_SP_SP

end module map
