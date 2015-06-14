!##############################################################################
!# Tutorial 004f: Datastructures - comparison of linked lists vs. maps
!##############################################################################

module tutorial004f

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use listInt
  use mapInt
  use random
  use statistics

  implicit none
  private

  public :: start_tutorial004f

contains

  ! ***************************************************************************

  subroutine start_tutorial004f

    ! Declare some variables
    type(t_random) :: rrandom
    type(t_timer) :: rtimer
    type(t_listInt) :: rlist
    type(it_listInt) :: rlistIterator
    type(t_mapInt) :: rmap
    type(it_mapInt) :: rmapIterator
    integer(I32), dimension(:), allocatable :: Irandom
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004f")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Initialise random number generator
    ! =================================

    call rng_initByClock(rrandom)

    ! =================================
    ! Generate 20000 random integer values
    ! =================================

    allocate(Irandom(20000))
    do i=1,size(Irandom)
      call rng_get(rrandom,Irandom(i))
    end do

    ! =================================
    ! Start time measurement for list
    ! =================================

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! =================================
    ! Create list
    ! =================================

    call list_create(rlist, size(Irandom))

    ! =================================
    ! Insert items at the end
    ! =================================

    do i=1,size(Irandom)
#ifdef USE_LARGEINT
      call list_push_back(rlist,int(Irandom(i), I64))
#else
      call list_push_back(rlist,Irandom(i))
#endif
    end do

    ! =================================
    ! Search for all items in list
    ! =================================

    do i=1,size(Irandom)
#ifdef USE_LARGEINT
      rlistIterator = list_find(rlist, int(Irandom(i), I64))
#else
      rlistIterator = list_find(rlist, Irandom(i))
#endif
    end do

    ! =================================
    ! Release list
    ! =================================

    call list_release(rlist)

    ! =================================
    ! Stop time measurement
    ! =================================

    call stat_stopTimer(rtimer)
    call output_line ("Total time for linked list: "&
        // trim(adjustl(sys_sd(rtimer%delapsedReal,10)))&
        // " seconds")

    ! =================================
    ! Start time measurement for map
    ! =================================

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! =================================
    ! Create map
    ! =================================

    call map_create(rmap, size(Irandom))

    ! =================================
    ! Insert items
    ! =================================

    do i=1,size(Irandom)
#ifdef USE_LARGEINT
      rmapIterator = map_insert(rmap,int(Irandom(i), I64))
#else
      rmapIterator = map_insert(rmap,Irandom(i))
#endif
    end do

    ! =================================
    ! Search for all items in map
    ! =================================

    do i=1,size(Irandom)
#ifdef USE_LARGEINT
      rmapIterator = map_find(rmap, int(Irandom(i), I64))
#else
      rmapIterator = map_find(rmap, Irandom(i))
#endif
    end do

    ! =================================
    ! Release map
    ! =================================

    call map_release(rmap)

    ! =================================
    ! Stop time measurement
    ! =================================

    call stat_stopTimer(rtimer)
    call output_line ("Total time for map        : "&
        // trim(adjustl(sys_sd(rtimer%delapsedReal,10)))&
        // " seconds")

    ! =================================
    ! Free memory
    ! =================================

    deallocate(Irandom)

  end subroutine

end module
