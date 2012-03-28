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
  use octreeDble, only : t_octreeDble, otree_release
  use octreeSngl, only : t_octreeSngl, otree_release
  use storage

  implicit none

  private
  public :: octree_init
  public :: octree_done
  public :: octree_getbase

  interface octree_getbase
    module procedure octree_getbase_dble
    module procedure octree_getbase_sngl
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
    type(t_octreeDble), pointer :: p_octreeDble => null()

    ! Pointer to single-valued octree implementations
    type(t_octreeSngl), pointer :: p_octreeSngl => null()

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
      allocate(roctree%p_octreeDble)

    case (OCTREE_SINGLE)
      allocate(roctree%p_octreeSngl)

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

    if (associated(roctree%p_octreeDble)) then
      call otree_release(roctree%p_octreeDble)
      deallocate(roctree%p_octreeDble)
    end if

    if (associated(roctree%p_octreeSngl)) then
      call otree_release(roctree%p_octreeSngl)
      deallocate(roctree%p_octreeSngl)
    end if

  end subroutine octree_done

  !************************************************************************

!<subroutine>

  subroutine octree_getbase_dble(roctree, p_roctree)

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
    type(t_octreeDble), pointer :: p_roctree
!</output>
!</subroutine>

    if (associated(roctree%p_octreeDble)) then
      p_roctree => roctree%p_octreeDble
    else
      nullify(p_roctree)
    end if

  end subroutine octree_getbase_dble

  !************************************************************************

!<subroutine>

  subroutine octree_getbase_sngl(roctree, p_roctree)

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
    type(t_octreeSngl), pointer :: p_roctree
!</output>
!</subroutine>

    if (associated(roctree%p_octreeSngl)) then
      p_roctree => roctree%p_octreeSngl
    else
      nullify(p_roctree)
    end if

  end subroutine octree_getbase_sngl

end module octree
