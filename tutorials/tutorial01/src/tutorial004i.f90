!##############################################################################
!# Tutorial 004i: Datastructures - quadtrees
!##############################################################################

module tutorial004i

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use quadtreebase
  use quadtreeDP

  implicit none
  private
  
  public :: start_tutorial004i

contains

  ! ***************************************************************************

  subroutine start_tutorial004i

    ! Declare some variables
    type(t_quadtreeDP) :: rquadtree
    integer :: i,iresult

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004i")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create quadtree with bounding box (0,1)x(0,2)
    ! which can store 20 vertices in 10 nodes before
    ! new memory is allocated
    ! =================================

    call qtree_create(rquadtree, 20, 10, 0.0_DP, 0.0_DP, 1.0_DP, 2.0_DP)

    ! =================================
    ! Insert new vertices into quadtree
    ! =================================

    iresult = qtree_insert(rquadtree, (/0.2_DP, 0.6_DP/))
    iresult = qtree_insert(rquadtree, (/0.2_DP, 0.5_DP/))
    iresult = qtree_insert(rquadtree, (/0.3_DP, 0.6_DP/))
    iresult = qtree_insert(rquadtree, (/0.5_DP, 0.6_DP/))
    iresult = qtree_insert(rquadtree, (/0.6_DP, 0.7_DP/))
    iresult = qtree_insert(rquadtree, (/0.8_DP, 2.0_DP/))
    iresult = qtree_insert(rquadtree, (/0.7_DP, 0.8_DP/))
    iresult = qtree_insert(rquadtree, (/0.4_DP, 1.6_DP/))

    select case(iresult)
    case (QTREE_FAILED)
      call output_line ("Failed to insert vertex into quadtree!")
    case (QTREE_FOUND)
      call output_line ("Vertex is already present in quadtree!")
    case (QTREE_INSERTED)
      call output_line ("Inserted vertex into quadtree!")
    end select

    iresult = qtree_insert(rquadtree, (/0.4_DP, 1.6_DP/))

    select case(iresult)
    case (QTREE_FAILED)
      call output_line ("Failed to insert vertex into quadtree!")
    case (QTREE_FOUND)
      call output_line ("Vertex is already present in quadtree!")
    case (QTREE_INSERTED)
      call output_line ("Inserted vertex into quadtree!")
    end select

    ! =================================
    ! Delete vertex from quadtree
    ! =================================

    iresult = qtree_delete(rquadtree, (/0.2_DP, 0.6_DP/))

    select case(iresult)
    case (QTREE_FAILED)
      call output_line ("Failed to delete vertex from quadtree!")
    case (QTREE_DELETED)
      call output_line ("Deleted vertex from quadtree!")
    end select

    iresult = qtree_delete(rquadtree, (/0.1_DP, 1.3_DP/))

    select case(iresult)
    case (QTREE_FAILED)
      call output_line ("Failed to delete vertex from quadtree!")
    case (QTREE_DELETED)
      call output_line ("Deleted vertex from quadtree!")
    end select

    ! =================================
    ! Find vertex in quadtree
    ! =================================

    iresult = qtree_find(rquadtree, (/0.501_DP, 0.601_DP/))

    select case(iresult)
    case (QTREE_NOT_FOUND)
      call output_line ("Unable to find vertex in quadtree!")
    case (QTREE_FOUND)
      call output_line ("Found vertex in quadtree!")
    end select

    ! =================================
    ! Find vertex in quadtree using user-defined isEqual
    ! function (with less restrictive proximity definition)
    ! =================================

    iresult = qtree_find(rquadtree, (/0.501_DP, 0.601_DP/),&
                         fcb_isEqual=my_isEqual)

    select case(iresult)
    case (QTREE_NOT_FOUND)
      call output_line ("Unable to find vertex in quadtree!")
    case (QTREE_FOUND)
      call output_line ("Found vertex in quadtree!")
    end select

    ! =================================
    ! Release quadtree
    ! =================================

    call qtree_release(rquadtree)

  end subroutine

  ! ***************************************************************************

  pure logical function my_isEqual(data1,data2) result(bisEqual)
    real(DP), dimension(2), intent(in) :: data1,data2

    bisEqual = (maxval(abs(data1-data2)) .lt. 1e-1)

  end function my_isEqual
  
end module
