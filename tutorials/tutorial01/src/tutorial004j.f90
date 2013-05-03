!##############################################################################
!# Tutorial 004j: Datastructures - octrees
!##############################################################################

module tutorial004j

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use octreebase
  use octreeSP

  implicit none
  private
  
  public :: start_tutorial004j

contains

  ! ***************************************************************************

  subroutine start_tutorial004j

    ! Declare some variables
    type(t_octreeSP) :: roctree
    integer :: i,iresult

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004j")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create octree with bounding box (0,1)x(0,2)
    ! which can store 20 vertices in 10 nodes before
    ! new memory is allocated
    ! =================================

    call otree_create(roctree, 20, 10, 0.0_SP, 0.0_SP, 0.0_SP, 1.0_SP, 2.0_SP, 1.4_SP)

    ! =================================
    ! Insert new vertices into octree
    ! =================================

    iresult = otree_insert(roctree, (/0.2_SP, 0.6_SP, 1.1_SP/))
    iresult = otree_insert(roctree, (/0.2_SP, 0.5_SP, 0.3_SP/))
    iresult = otree_insert(roctree, (/0.3_SP, 0.6_SP, 0.3_SP/))
    iresult = otree_insert(roctree, (/0.5_SP, 0.6_SP, 0.7_SP/))
    iresult = otree_insert(roctree, (/0.6_SP, 0.7_SP, 1.0_SP/))
    iresult = otree_insert(roctree, (/0.8_SP, 2.0_SP, 0.8_SP/))
    iresult = otree_insert(roctree, (/0.7_SP, 0.8_SP, 0.7_SP/))
    iresult = otree_insert(roctree, (/0.4_SP, 1.6_SP, 0.4_SP/))

    select case(iresult)
    case (OTREE_FAILED)
      call output_line ("Failed to insert vertex into octree!")
    case (OTREE_FOUND)
      call output_line ("Vertex is already present in octree!")
    case (OTREE_INSERTED)
      call output_line ("Inserted vertex into octree!")
    end select

    iresult = otree_insert(roctree, (/0.4_SP, 1.6_SP, 0.4_SP/))

    select case(iresult)
    case (OTREE_FAILED)
      call output_line ("Failed to insert vertex into octree!")
    case (OTREE_FOUND)
      call output_line ("Vertex is already present in octree!")
    case (OTREE_INSERTED)
      call output_line ("Inserted vertex into octree!")
    end select

    ! =================================
    ! Delete vertex from octree
    ! =================================

    iresult = otree_delete(roctree, (/0.2_SP, 0.6_SP, 1.1_SP/))

    select case(iresult)
    case (OTREE_FAILED)
      call output_line ("Failed to delete vertex from octree!")
    case (OTREE_DELETED)
      call output_line ("Deleted vertex from octree!")
    end select

    iresult = otree_delete(roctree, (/0.1_SP, 1.3_SP, 1.1_SP/))

    select case(iresult)
    case (OTREE_FAILED)
      call output_line ("Failed to delete vertex from octree!")
    case (OTREE_DELETED)
      call output_line ("Deleted vertex from octree!")
    end select

    ! =================================
    ! Find vertex in octree
    ! =================================

    iresult = otree_find(roctree, (/0.501_SP, 0.6_SP, 0.7_SP/))

    select case(iresult)
    case (OTREE_NOT_FOUND)
      call output_line ("Unable to find vertex in octree!")
    case (OTREE_FOUND)
      call output_line ("Found vertex in octree!")
    end select

    ! =================================
    ! Find vertex in octree using user-defined isEqual
    ! function (with less restrictive proximity definition)
    ! =================================

    iresult = otree_find(roctree, (/0.501_SP, 0.6_SP, 0.7_SP/),&
                         fcb_isEqual=my_isEqual)

    select case(iresult)
    case (OTREE_NOT_FOUND)
      call output_line ("Unable to find vertex in octree!")
    case (OTREE_FOUND)
      call output_line ("Found vertex in octree!")
    end select

    ! =================================
    ! Release octree
    ! =================================

    call otree_release(roctree)

  end subroutine

  ! ***************************************************************************

  pure logical function my_isEqual(data1,data2) result(bisEqual)
    real(SP), dimension(3), intent(in) :: data1,data2

    bisEqual = (maxval(abs(data1-data2)) .lt. 1e-1)

  end function my_isEqual
  
end module
