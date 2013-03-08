!##############################################################################
!# ****************************************************************************
!# <name> octree </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a meta octree, that is, this data structure provides
!# pointers to all standard octree implementations available.
!#
!# The following routine are available:
!#
!# 1.) octree_init
!#     -> Initialises a meta octree structure
!#
!# 2.) octree_done
!#     -> Finalises a meta octree structure
!#
!# 3.) octree_getbase
!#     -> Returns pointer to a concrete octree
!#
!# </purpose>
!##############################################################################

module octree

!$use omp_lib
  use fsystem
  use genoutput
  use octreeDP, only : t_octreeDP, otree_release
  use octreeSP, only : t_octreeSP, otree_release
  use storage

  implicit none

  private
  public :: t_octree
  public :: octree_init
  public :: octree_done
  public :: octree_getbase

  interface octree_getbase
    module procedure octree_getbase_DP
    module procedure octree_getbase_SP
  end interface

!<constants>
!<constantblock description="Global flags for octree implementations">

  ! octree for double data
  integer, parameter, public :: OCTREE_DOUBLE = ST_DOUBLE

  ! octree for single data
  integer, parameter, public :: OCTREE_SINGLE = ST_SINGLE

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! Meta octree structure
  type t_octree
    private

    ! Pointer to double-valued octree implementations
    type(t_octreeDP), pointer :: p_octreeDP => null()

    ! Pointer to single-valued octree implementations
    type(t_octreeSP), pointer :: p_octreeSP => null()

  end type t_octree

!</typeblock>
!</types>

contains

  !************************************************************************

!<subroutine>

  subroutine octree_init(roctree, coctreeType)

!<description>
    ! Initialises a meta octree structure
!</description>

!<input>
    ! octree type
    integer, intent(in) :: coctreeType
!</input>

!<output>
    ! Meta octree
    type(t_octree), intent(out) :: roctree
!</output>
!</subroutine>

    select case (coctreeType)    
    case (OCTREE_DOUBLE)
      allocate(roctree%p_octreeDP)

    case (OCTREE_SINGLE)
      allocate(roctree%p_octreeSP)

    case DEFAULT
      call output_line('Invalid octree type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'octree_init')
      call sys_halt()
    end select

  end subroutine octree_init

  !************************************************************************

!<subroutine>

  subroutine octree_done(roctree)

!<description>
    ! Finalises a meta octree structure
!</description>

!<inputoutput>
    ! Meta octree
    type(t_octree), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

    if (associated(roctree%p_octreeDP)) then
      call otree_release(roctree%p_octreeDP)
      deallocate(roctree%p_octreeDP)
    end if

    if (associated(roctree%p_octreeSP)) then
      call otree_release(roctree%p_octreeSP)
      deallocate(roctree%p_octreeSP)
    end if

  end subroutine octree_done

  !************************************************************************

!<subroutine>

  subroutine octree_getbase_DP(roctree, p_roctree)

!<description>
    ! Returns a pointer to the double-valued octree implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta octree
    type(t_octree), intent(in) :: roctree
!</input>

!<output>
    ! Pointer to the octree implementation
    type(t_octreeDP), pointer :: p_roctree
!</output>
!</subroutine>

    if (associated(roctree%p_octreeDP)) then
      p_roctree => roctree%p_octreeDP
    else
      nullify(p_roctree)
    end if

  end subroutine octree_getbase_DP

  !************************************************************************

!<subroutine>

  subroutine octree_getbase_SP(roctree, p_roctree)

!<description>
    ! Returns a pointer to the single-valued octree implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta octree
    type(t_octree), intent(in) :: roctree
!</input>

!<output>
    ! Pointer to the octree implementation
    type(t_octreeSP), pointer :: p_roctree
!</output>
!</subroutine>

    if (associated(roctree%p_octreeSP)) then
      p_roctree => roctree%p_octreeSP
    else
      nullify(p_roctree)
    end if

  end subroutine octree_getbase_SP

end module octree
