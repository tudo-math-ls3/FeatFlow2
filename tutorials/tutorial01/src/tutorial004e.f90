!##############################################################################
!# Tutorial 004e: Datastructures - maps
!##############################################################################

module tutorial004e

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use mapInt

  implicit none
  private
  
  public :: start_tutorial004e

contains

  ! ***************************************************************************

  subroutine start_tutorial004e

    ! Declare some variables
    type(t_mapInt) :: rmap
    type(it_mapInt) :: riterator
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004e")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create map with initial space for 20
    ! =================================

    call map_create(rmap, 20)

    ! =================================
    ! Insert new items into the map
    ! =================================

    riterator = map_insert(rmap, 5)
    riterator = map_insert(rmap, 23)
    riterator = map_insert(rmap, 2)
    riterator = map_insert(rmap, 76)
    riterator = map_insert(rmap, 9)
    riterator = map_insert(rmap, 9)
    riterator = map_insert(rmap, -2)
    riterator = map_insert(rmap, -52)
    riterator = map_insert(rmap, -4)
    riterator = map_insert(rmap, -52)

    ! =================================
    ! Check size and content of the arraylist
    ! =================================
    i = map_size(rmap)
    call output_line ("Size of the map: " // trim(sys_siL(i,10)))
    i = map_max_size(rmap)
    call output_line ("Maximum size of the map (before internal data is reallocated): "&
        // trim(sys_siL(i,10)))

    if (map_empty(rmap)) then
      call output_line ("Map is empty!")
    else 
      call output_line ("Map is not empty!")
    end if

    ! =================================
    ! Perform forward iteration through map items
    ! =================================

    riterator = map_begin(rmap)
    do while (riterator .ne. map_end(rmap))
      i = map_get(rmap, riterator)
      call output_line (" " // trim(sys_siL(i,10)), bnolinebreak=.true.)
      call map_next(riterator)     
    end do
    call output_lbrk()
    
    ! =================================
    ! Search for item in map
    ! =================================

    riterator = map_find(rmap, 2)
    if (map_isNull(riterator)) then
      call output_line ("Map does not contain item 2!")
    else
      call output_line ("Map contains item 2!")
    end if

    ! =================================
    ! Erase item from map
    ! =================================

    call map_erase(rmap, riterator)

    ! =================================
    ! Search for item in map
    ! =================================

    riterator = map_find(rmap, 2)
    if (map_isNull(riterator)) then
      call output_line ("Map does not contain item 2!")
    else
      call output_line ("Map contains item 2!")
    end if

    ! =================================
    ! Release map
    ! =================================

    call map_release(rmap)

  end subroutine

end module
