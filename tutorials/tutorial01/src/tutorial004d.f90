!##############################################################################
!# Tutorial 004d: Datastructures - linked lists advanced features
!##############################################################################

module tutorial004d

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use listInt

  implicit none
  private
  
  public :: start_tutorial004d

contains

  ! ***************************************************************************

  subroutine start_tutorial004d

    ! Declare some variables
    type(t_listInt) :: rlist1,rlist2
    type(it_listInt) :: riterator
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004d")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create two linked list with initial space for 10 items
    ! =================================

    call list_create(rlist1, 10)
    call list_create(rlist2, 10)

    ! =================================
    ! Insert new items at the end of the list #1
    ! =================================

    call list_push_back(rlist1,  5)
    call list_push_back(rlist1,  23)
    call list_push_back(rlist1,  2)
    call list_push_back(rlist1, -2)
    call list_push_back(rlist1,  76)
    call list_push_back(rlist1,  9)

    ! =================================
    ! Insert new items at the beginning of the list #2
    ! =================================

    call list_push_front(rlist2,  9)
    call list_push_front(rlist2, -2)
    call list_push_front(rlist2, -52)
    call list_push_front(rlist2, -4)
    call list_push_front(rlist2,  2)
    call list_push_front(rlist2,  5)

    ! =================================
    ! Sort both lists
    ! =================================

    call list_sort(rlist1)
    call list_sort(rlist2)

    ! =================================
    ! Merge list #1 into list #2
    ! =================================
    
    call list_merge(rlist2, rlist1)

    ! =================================
    ! Perform forward iteration through list #2
    ! =================================

    riterator = list_begin(rlist2)
    do while (riterator .ne. list_end(rlist2))
      i = list_get(rlist2, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call list_next(riterator)     
    end do
    call output_lbrk()

    ! =================================
    ! Remove duplicate items from sorted list
    ! =================================

    call list_unique(rlist2)

    ! =================================
    ! Perform forward iteration through list #2
    ! =================================

    riterator = list_begin(rlist2)
    do while (riterator .ne. list_end(rlist2))
      i = list_get(rlist2, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call list_next(riterator)     
    end do
    call output_lbrk()

    ! =================================
    ! Release linked lists
    ! =================================

    call list_release(rlist1)
    call list_release(rlist2)

  end subroutine

end module
