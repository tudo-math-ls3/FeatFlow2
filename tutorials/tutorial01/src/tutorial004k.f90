!##############################################################################
!# Tutorial 004k: Datastructures - objects
!##############################################################################

module tutorial004k

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use listDP
  use mapInt
  use quadtreeSP

  implicit none
  private
  
  public :: start_tutorial004k

contains

  ! ***************************************************************************

  subroutine start_tutorial004k

    ! Declare some variables
    type(t_genericObject) :: robject
    type(t_listDP) :: rlist
    type(t_mapInt) :: rmap
    type(it_mapInt) :: rmapIterator
    type(t_quadtreeSP) :: rquadtree
    integer :: iresult

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004k")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create list, map, and quadtree
    ! =================================

    call list_create(rlist, 10)
    call map_create(rmap, 10)
    call qtree_create(rquadtree, 20, 10, 0.0_SP, 0.0_SP, 1.0_SP, 2.0_SP)

    ! =================================
    ! Insert new items into the list
    ! =================================

    call list_push_back(rlist,  1.2_DP)
    call list_push_back(rlist,  3.3_DP)
    call list_push_back(rlist, 34.0_DP)
    call list_push_back(rlist,  4.5_DP)
    call list_push_back(rlist,  9.4_DP)
    call list_push_back(rlist, 14.5_DP)
    call list_push_back(rlist, -4.2_DP)
    call list_push_back(rlist,  8.7_DP)

    ! =================================
    ! Insert new items into the map
    ! =================================
    
    rmapIterator = map_insert(rmap, 5)
    rmapIterator = map_insert(rmap, 23)
    rmapIterator = map_insert(rmap, 2)
    rmapIterator = map_insert(rmap, 76)
    rmapIterator = map_insert(rmap, 9)
    rmapIterator = map_insert(rmap, 9)
    rmapIterator = map_insert(rmap, -2)
    rmapIterator = map_insert(rmap, -52)
    rmapIterator = map_insert(rmap, -4)
    rmapIterator = map_insert(rmap, -52)
    
    ! =================================
    ! Insert new vertices into quadtree
    ! =================================

    iresult = qtree_insert(rquadtree, (/0.2_SP, 0.6_SP/))
    iresult = qtree_insert(rquadtree, (/0.2_SP, 0.5_SP/))
    iresult = qtree_insert(rquadtree, (/0.3_SP, 0.6_SP/))
    iresult = qtree_insert(rquadtree, (/0.5_SP, 0.6_SP/))
    iresult = qtree_insert(rquadtree, (/0.6_SP, 0.7_SP/))
    iresult = qtree_insert(rquadtree, (/0.8_SP, 2.0_SP/))
    iresult = qtree_insert(rquadtree, (/0.7_SP, 0.8_SP/))
    iresult = qtree_insert(rquadtree, (/0.4_SP, 1.6_SP/))
    
    ! =================================
    ! Cast list, map, and quadtree to generic object
    ! and pass it as argument to print routine
    ! =================================    

    call list_cast(rlist, robject)
    call output_line("Print content of linked list")
    call print_genericObject(robject,1)
    call output_lbrk()

    call map_cast(rmap, robject)
    call output_line("Print content of map")
    call print_genericObject(robject,2)
    call output_lbrk()

    call qtree_cast(rquadtree, robject)
    call output_line("Print content of quadtree")
    call print_genericObject(robject,3)
    call output_lbrk()
    
    ! =================================
    ! Release data structures
    ! =================================

    call list_release(rlist)
    call map_release(rmap)
    call qtree_release(rquadtree)

  end subroutine

  ! ***************************************************************************

  subroutine print_genericObject(robject,ctype)
    type(t_genericObject), intent(in) :: robject
    integer, intent(in) :: ctype
    character(LEN=SYS_STRLEN) :: spostdir

    type(t_listDP), pointer :: p_rlist
    type(t_mapInt), pointer :: p_rmap
    type(t_quadtreeSP), pointer :: p_rquadtree

    select case(ctype)
    case(1)
      call list_uncast(robject, p_rlist)
      call list_print(p_rlist)
    case(2)
      call map_uncast(robject, p_rmap)
      call map_print(p_rmap)
    case(3)
      call qtree_uncast(robject, p_rquadtree)
      if (sys_getenv_string("POSTDIR",spostdir)) then
        call qtree_print(p_rquadtree, trim(spostdir)//"/tutorial004k_quadtree.mat")
      else
        call qtree_print(p_rquadtree, "post/tutorial004k_quadtree.mat")
      end if
      call output_line("The quadtree can be visualized by running the Matlab script")
      call output_line("> plotquadtree('post/tutorial004k_quadtree.mat')")
      call output_line("Note that you need to have Featflow2/matlab in the path.")
      
    case default
      call output_line ("This given object type is not supported!")
      call sys_halt()
    end select
    
  end subroutine
  
end module
