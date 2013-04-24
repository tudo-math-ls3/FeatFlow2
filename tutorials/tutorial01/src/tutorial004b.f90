!##############################################################################
!# Tutorial 004b: Datastructures - linked lists
!##############################################################################

module tutorial004b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use listInt
  use sort

  implicit none
  private
  
  public :: start_tutorial004b

contains

  ! ***************************************************************************

  subroutine start_tutorial004b

    ! Declare some variables
    type(t_listInt) :: rlist
    type(it_listInt) :: riterator
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004b")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create linked list
    ! =================================

    call list_create(rlist, 10)

    ! =================================
    ! Insert new items at the end of the list
    ! =================================

    call list_push_back(rlist, 5)
    call list_push_back(rlist, 23)
    call list_push_back(rlist, 2)
    call list_push_back(rlist, 76)

    ! =================================
    ! Insert new items at the beginning of the list
    ! =================================

    call list_push_front(rlist,  9)
    call list_push_front(rlist, -2)
    call list_push_front(rlist, -52)
    call list_push_front(rlist, -4)

    ! =================================
    ! Check size and content of the list
    ! =================================
    i = list_size(rlist)
    call output_line ("Size of the list: " // trim(sys_siL(i,10)))
    i = list_max_size(rlist)
    call output_line ("Maximum size of the list (before internal data is reallocated): "&
        // trim(sys_siL(i,10)))
    if (list_empty(rlist)) then
      call output_line ("List is empty!")
    else 
      call output_line ("List is not empty!")
    end if

    ! =================================
    ! Perform forward iteration through list items
    ! =================================

    riterator = list_begin(rlist)
    do while (riterator .ne. list_end(rlist))
      i = list_get(rlist, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call list_next(riterator)     
    end do
    call output_lbrk()
    
    ! =================================
    ! Perform backward iteration through list items
    ! =================================

    riterator = list_rbegin(rlist)
    do while (riterator .ne. list_rend(rlist))
      i = list_get(rlist, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call list_next(riterator)     
    end do
    call output_lbrk()

    ! =================================
    ! Reverse list items
    ! =================================

    call list_reverse(rlist)

    ! =================================
    ! Perform forward iteration through list items
    ! =================================

    riterator = list_begin(rlist)
    do while (riterator .ne. list_end(rlist))
      i = list_get(rlist, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call list_next(riterator)     
    end do
    call output_lbrk()

    ! =================================
    ! Sort list items
    ! =================================

    call list_sort(rlist)

    ! =================================
    ! Perform forward iteration through list items
    ! =================================

    riterator = list_begin(rlist)
    do while (riterator .ne. list_end(rlist))
      i = list_get(rlist, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call list_next(riterator)     
    end do
    call output_lbrk()

    ! =================================
    ! Remove first and last item from list
    ! =================================
    call list_pop_front(rlist)
    call list_pop_back(rlist)

    ! =================================
    ! Perform forward iteration through list items
    ! =================================
    
    riterator = list_begin(rlist)
    do while (riterator .ne. list_end(rlist))
      i = list_get(rlist, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call list_next(riterator)     
    end do
    call output_lbrk()

    ! =================================
    ! Find item by key value and erase it
    ! =================================
    riterator = list_find(rlist, 9)
    if (.not.list_isNull(riterator)) then
      riterator = list_erase(rlist, riterator)
    end if

    ! =================================
    ! Insert item it position
    ! =================================
    
    riterator = list_insert(rlist, riterator, 47)

    ! =================================
    ! Perform forward iteration through list items
    ! =================================
    
    riterator = list_begin(rlist)
    do while (riterator .ne. list_end(rlist))
      i = list_get(rlist, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call list_next(riterator)     
    end do
    call output_lbrk()

    ! =================================
    ! Release linked list
    ! =================================

    call list_release(rlist)

    stop

  end subroutine

end module
