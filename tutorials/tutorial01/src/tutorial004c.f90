!##############################################################################
!# Tutorial 004c: Datastructures - arraylists
!##############################################################################

module tutorial004c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use arraylistInt

  implicit none
  private
  
  public :: start_tutorial004c

contains

  ! ***************************************************************************

  subroutine start_tutorial004c

    ! Declare some variables
    type(t_arraylistInt) :: rarraylist
    type(it_arraylistInt) :: riterator
    integer :: i,j

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004c")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create arraylist with 5 tables and initial space for 10
    ! =================================

    call alst_create(rarraylist, 5, 50)

    ! =================================
    ! Insert new items at the end of each list in tables [1,2,3,4]
    ! =================================

    call alst_push_back(rarraylist, 1, 5)
    call alst_push_back(rarraylist, 2, 23)
    call alst_push_back(rarraylist, 3, 2)
    call alst_push_back(rarraylist, 4, 76)
    call alst_push_back(rarraylist, 1, 13)
    call alst_push_back(rarraylist, 3, 1)
    call alst_push_back(rarraylist, 2, 44)
    call alst_push_back(rarraylist, 1, 73)
    call alst_push_back(rarraylist, 4, 11)
    call alst_push_back(rarraylist, 4, 53)
    call alst_push_back(rarraylist, 2, 2)
    call alst_push_back(rarraylist, 1, 9)

    ! =================================
    ! Insert new items at the beginning of each list in tables [1,2,3,4]
    ! =================================

    call alst_push_front(rarraylist, 2,  9)
    call alst_push_front(rarraylist, 4, -2)
    call alst_push_front(rarraylist, 2, -52)
    call alst_push_front(rarraylist, 4, -4)
    call alst_push_front(rarraylist, 1, -9)
    call alst_push_front(rarraylist, 3, -22)
    call alst_push_front(rarraylist, 3, -12)
    call alst_push_front(rarraylist, 4, -7)

    ! =================================
    ! Check size and content of the arraylist
    ! =================================
    i = alst_ntable(rarraylist)
    call output_line ("Number of tables in the arraylist: " // trim(sys_siL(i,10)))
    i = alst_size(rarraylist)
    call output_line ("Size of the arraylist: " // trim(sys_siL(i,10)))
    i = alst_max_size(rarraylist)
    call output_line ("Maximum size of the arraylist (before internal data is reallocated): "&
        // trim(sys_siL(i,10)))

    if (alst_empty(rarraylist)) then
      call output_line ("Arraylist is empty!")
    else 
      call output_line ("Arraylist is not empty!")
    end if

    do i=1,alst_ntable(rarraylist)
      if (alst_empty(rarraylist,i)) then
        call output_line ("List " // trim(sys_siL(i,10)) // " is empty!")
      else 
        call output_line ("List " // trim(sys_siL(i,10)) // " is not empty!")
      end if
    end do

    ! =================================
    ! Perform forward iteration through all(!) list items
    ! =================================

    do i=1,alst_ntable(rarraylist)
      call output_line ("Table " // trim(sys_siL(i,10)) // " : ", bnolinebreak=.true.)
      riterator = alst_begin(rarraylist,i)
      do while (riterator .ne. alst_end(rarraylist,i))
        j = alst_get(rarraylist, riterator)
        call output_line (" " // trim(sys_siL(j,10)), bnolinebreak=.true.)
        call alst_next(riterator)     
      end do
      call output_lbrk()
    end do
    
    ! =================================
    ! Perform backward iteration through all(!) list items
    ! =================================

    do i=1,alst_ntable(rarraylist)
      call output_line ("Table " // trim(sys_siL(i,10)) // " : ", bnolinebreak=.true.)
      riterator = alst_rbegin(rarraylist,i)
      do while (riterator .ne. alst_rend(rarraylist,i))
        j = alst_get(rarraylist, riterator)
        call output_line (" " // trim(sys_siL(j,10)), bnolinebreak=.true.)
        call alst_next(riterator)     
      end do
      call output_lbrk()
    end do

    ! =================================
    ! Reverse all(!) list items in all tables
    ! =================================

    call alst_reverse(rarraylist)

    ! =================================
    ! Perform forward iteration through all(!) list items
    ! =================================

    do i=1,alst_ntable(rarraylist)
      call output_line ("Table " // trim(sys_siL(i,10)) // " : ", bnolinebreak=.true.)
      riterator = alst_begin(rarraylist,i)
      do while (riterator .ne. alst_end(rarraylist,i))
        j = alst_get(rarraylist, riterator)
        call output_line (" " // trim(sys_siL(j,10)), bnolinebreak=.true.)
        call alst_next(riterator)     
      end do
      call output_lbrk()
    end do

    ! =================================
    ! Sort all(!) list items
    ! =================================

    call alst_sort(rarraylist)

    ! =================================
    ! Reverse list items in table 2
    ! =================================

    call alst_reverse(rarraylist,2)

    ! =================================
    ! Perform forward iteration through all(!) list items
    ! =================================

    do i=1,alst_ntable(rarraylist)
      call output_line ("Table " // trim(sys_siL(i,10)) // " : ", bnolinebreak=.true.)
      riterator = alst_begin(rarraylist,i)
      do while (riterator .ne. alst_end(rarraylist,i))
        j = alst_get(rarraylist, riterator)
        call output_line (" " // trim(sys_siL(j,10)), bnolinebreak=.true.)
        call alst_next(riterator)     
      end do
      call output_lbrk()
    end do

    ! =================================
    ! Clear second table and remove last item from third table
    ! =================================
    call alst_clear(rarraylist,2)
    call alst_pop_back(rarraylist,3)

    ! =================================
    ! Perform forward iteration through all(!) list items
    ! =================================
    
    do i=1,alst_ntable(rarraylist)
      call output_line ("Table " // trim(sys_siL(i,10)) // " : ", bnolinebreak=.true.)
      riterator = alst_begin(rarraylist,i)
      do while (riterator .ne. alst_end(rarraylist,i))
        j = alst_get(rarraylist, riterator)
        call output_line (" " // trim(sys_siL(j,10)), bnolinebreak=.true.)
        call alst_next(riterator)     
      end do
      call output_lbrk()
    end do

    ! =================================
    ! Find item by key value and erase it
    ! =================================
    riterator = alst_find(rarraylist, 3, -12)
    if (.not.alst_isNull(riterator)) then
      riterator = alst_erase(rarraylist, riterator)
    end if

    ! =================================
    ! Insert item at position
    ! =================================
    
    riterator = alst_insert(rarraylist, riterator, 47)

    ! =================================
    ! Perform forward iteration through all(!) list items
    ! =================================
    
    do i=1,alst_ntable(rarraylist)
      call output_line ("Table " // trim(sys_siL(i,10)) // " : ", bnolinebreak=.true.)
      riterator = alst_begin(rarraylist,i)
      do while (riterator .ne. alst_end(rarraylist,i))
        j = alst_get(rarraylist, riterator)
        call output_line (" " // trim(sys_siL(j,10)), bnolinebreak=.true.)
        call alst_next(riterator)     
      end do
      call output_lbrk()
    end do

    ! =================================
    ! Release table 3
    ! =================================

    call alst_release(rarraylist, 3)

    ! =================================
    ! Create new table 7 with item 42
    ! =================================

    call alst_create(rarraylist, 7)
    call alst_push_back(rarraylist, 7, 42)

    ! =================================
    ! Perform forward iteration through all(!) list items
    ! =================================
    
    do i=1,alst_ntable(rarraylist)
      call output_line ("Table " // trim(sys_siL(i,10)) // " : ", bnolinebreak=.true.)
      riterator = alst_begin(rarraylist,i)
      do while (riterator .ne. alst_end(rarraylist,i))
        j = alst_get(rarraylist, riterator)
        call output_line (" " // trim(sys_siL(j,10)), bnolinebreak=.true.)
        call alst_next(riterator)     
      end do
      call output_lbrk()
    end do

    ! =================================
    ! Release arraylist
    ! =================================

    call alst_release(rarraylist)

  end subroutine

end module
