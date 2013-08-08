!##############################################################################
!# Tutorial 004m: Datastructures - sorting of arrays
!##############################################################################

module tutorial004m

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use numbersets

  implicit none
  private
  
  public :: start_tutorial004m

contains

  ! ***************************************************************************

  subroutine start_tutorial004m

    ! Declare some variables
    type(t_directAccessIntSet) :: rset
    integer, dimension(50) :: Ilist
    integer :: nelementsInList,i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004m")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a set, put some numbers in
    ! it and print
    ! =================================
    
    ! Initialise a set
    call nsets_initDASet (rset,1040)
    
    ! Add some numbers
    call nsets_addToDASet (rset,1023)
    call nsets_addToDASet (rset,1024)
    call nsets_addToDASet (rset,1025)
    call nsets_addToDASet (rset,5)
    call nsets_addToDASet (rset,700)
    call nsets_addToDASet (rset,1002)
    call nsets_addToDASet (rset,1000)
    call nsets_addToDASet (rset,1000)
    call nsets_addToDASet (rset,1000)
    call nsets_addToDASet (rset,1001)
    call nsets_addToDASet (rset,1003)
    
    ! Remove some numbers
    call nsets_removeFromDASet (rset,1002)
    call nsets_removeFromDASet (rset,1025)
    
    ! Create a list with the elements in the set
    call nsets_getElements (rset,Ilist,nelementsInList)
    
    ! Print the list. The 1000 should be printed only once.
    call output_lbrk()
    do i=1,nelementsInList
      call output_line (sys_siL(Ilist(i),10))
    end do
    
    ! =================================
    ! Check if the list contains some
    ! elements
    ! =================================
    call output_lbrk()
    
    call output_line ("1000 in list: "//sys_siL( nsets_DASetContains(rset,1000) ,10))
    call output_line (" 800 in list: "//sys_siL( nsets_DASetContains(rset,800)  ,10))
    call output_line ("1023 in list: "//sys_siL( nsets_DASetContains(rset,1023) ,10))
    call output_line ("1024 in list: "//sys_siL( nsets_DASetContains(rset,1024) ,10))
    call output_line ("1025 in list: "//sys_siL( nsets_DASetContains(rset,1025) ,10))
    call output_line ("1030 in list: "//sys_siL( nsets_DASetContains(rset,1030) ,10))
    call output_line ("1035 in list: "//sys_siL( nsets_DASetContains(rset,1035) ,10))
    
    ! =================================
    ! Clear the list, add some other 
    ! elements
    ! =================================
    
    call nsets_clearDASet (rset)
    
    call nsets_putElements (rset, (/ 10, 11, 12, 50, 60, 100 /) )
    
    call output_lbrk()
    call output_line ("Elements in new list: " // sys_siL( nsets_nelementsInSet(rset), 10) )

    ! =================================
    ! Clean up the list
    ! =================================
    
    call nsets_doneDASet (rset)

  end subroutine

end module
