!##############################################################################
!# ****************************************************************************
!# <name> quadtree </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a meta quadtree, that is, this data structure provides
!# pointers to all standard quadtree implementations available.
!#
!# The following routine are available:
!#
!# 1.) quadtree_init
!#     -> Initialises a meta quadtree structure
!#
!# 2.) quadtree_done
!#     -> Finalises a meta quadtree structure
!#
!# 3.) quadtree_getbase
!#     -> Returns pointer to a concrete quadtree
!#
!# </purpose>
!##############################################################################

module quadtree

!$use omp_lib
  use fsystem
  use genoutput
  use quadtreeDP, only : t_quadtreeDP,qtree_release
  use quadtreeSP, only : t_quadtreeSP,qtree_release
  use storage

  implicit none

  private
  public :: t_quadtree
  public :: quadtree_init
  public :: quadtree_done
  public :: quadtree_getbase

  interface quadtree_getbase
    module procedure quadtree_getbase_DP
    module procedure quadtree_getbase_SP
  end interface

!<constants>
!<constantblock description="Global flags for quadtree implementations">

  ! quadtree for double data
  integer, parameter, public :: QUADTREE_DOUBLE = ST_DOUBLE

  ! quadtree for single data
  integer, parameter, public :: QUADTREE_SINGLE = ST_SINGLE

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! Meta quadtree structure
  type t_quadtree
    private

    ! Pointer to double-valued quadtree implementations
    type(t_quadtreeDP), pointer :: p_quadtreeDP => null()

    ! Pointer to single-valued quadtree implementations
    type(t_quadtreeSP), pointer :: p_quadtreeSP => null()

  end type t_quadtree

!</typeblock>
!</types>

contains

  !************************************************************************

!<subroutine>

  subroutine quadtree_init(rquadtree, cquadtreeType)

!<description>
    ! Initialises a meta quadtree structure
!</description>

!<input>
    ! quadtree type
    integer, intent(in) :: cquadtreeType
!</input>

!<output>
    ! Meta quadtree
    type(t_quadtree), intent(out) :: rquadtree
!</output>
!</subroutine>

    select case (cquadtreeType)    
    case (QUADTREE_DOUBLE)
      allocate(rquadtree%p_quadtreeDP)

    case (QUADTREE_SINGLE)
      allocate(rquadtree%p_quadtreeSP)

    case DEFAULT
      call output_line('Invalid quadtree type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'quadtree_init')
      call sys_halt()
    end select

  end subroutine quadtree_init

  !************************************************************************

!<subroutine>

  subroutine quadtree_done(rquadtree)

!<description>
    ! Finalises a meta quadtree structure
!</description>

!<inputoutput>
    ! Meta quadtree
    type(t_quadtree), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

    if (associated(rquadtree%p_quadtreeDP)) then
      call qtree_release(rquadtree%p_quadtreeDP)
      deallocate(rquadtree%p_quadtreeDP)
    end if

    if (associated(rquadtree%p_quadtreeSP)) then
      call qtree_release(rquadtree%p_quadtreeSP)
      deallocate(rquadtree%p_quadtreeSP)
    end if

  end subroutine quadtree_done

  !************************************************************************

!<subroutine>

  subroutine quadtree_getbase_DP(rquadtree, p_rquadtree)

!<description>
    ! Returns a pointer to the double-valued quadtree implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta quadtree
    type(t_quadtree), intent(in) :: rquadtree
!</input>

!<output>
    ! Pointer to the quadtree implementation
    type(t_quadtreeDP), pointer :: p_rquadtree
!</output>
!</subroutine>

    if (associated(rquadtree%p_quadtreeDP)) then
      p_rquadtree => rquadtree%p_quadtreeDP
    else
      nullify(p_rquadtree)
    end if

  end subroutine quadtree_getbase_DP

  !************************************************************************

!<subroutine>

  subroutine quadtree_getbase_SP(rquadtree, p_rquadtree)

!<description>
    ! Returns a pointer to the single-valued quadtree implementation
    ! without auxiliary data
!</description>

!<input>
    ! Meta quadtree
    type(t_quadtree), intent(in) :: rquadtree
!</input>

!<output>
    ! Pointer to the quadtree implementation
    type(t_quadtreeSP), pointer :: p_rquadtree
!</output>
!</subroutine>

    if (associated(rquadtree%p_quadtreeSP)) then
      p_rquadtree => rquadtree%p_quadtreeSP
    else
      nullify(p_rquadtree)
    end if

  end subroutine quadtree_getbase_SP

end module quadtree
