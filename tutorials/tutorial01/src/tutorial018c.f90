!##############################################################################
!# Tutorial 018c: Basic collection handling
!##############################################################################

module tutorial018c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use collection

  implicit none
  private
  
  public :: start_tutorial018c

contains

  ! ***************************************************************************

  subroutine printarray (rcollection)
  type(t_collection), intent(inout) :: rcollection
  
    integer, dimension(:), pointer :: p_Iarray
    integer :: i,itype
    logical :: bexists
    
    ! Does "MYFAIL" exist?
    bexists = collct_queryvalue (rcollection, "MYFAIL") .gt. 0
    call output_line ("MYFAIL   exists = "//sys_sl(bexists))
    
    ! Does "MYIARRAY" exist?
    bexists = collct_queryvalue (rcollection, "MYIARRAY") .gt. 0
    call output_line ("MYIARRAY exists = "//sys_sl(bexists))
    
    if (bexists) then
      ! Print the type
      itype = collct_gettype (rcollection, "MYIARRAY")
      call output_line ("Type            = " // sys_si(itype,5))

      ! Get a pointer to the data from the collection
      call collct_getvalue_iarrp (rcollection, "MYIARRAY", p_Iarray)
      
      ! Print the content
      do i=1,ubound(p_Iarray,1)
        call output_line (sys_si(p_Iarray(i),5))
      end do
    end if
  
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial018c

    ! Declare some variables.
    type(t_collection) :: rcollection
    integer, dimension(5), target :: Iarray
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 018c")
    call output_separator (OU_SEP_MINUS)

    ! For storing named items in the collection, the collection must be
    ! initialised.
    !
    ! Initialise the collection
    call collct_init (rcollection)
    
    ! Initialise an array
    Iarray(:) = (/ 5,4,3,2,1 /)
    
    ! Save a pointer to Iarray in the collection, named "MYIARRAY"
    call collct_setvalue_iarrp (rcollection, "MYIARRAY", Iarray, .true.)
    
    ! Pass the collection to the subroutine
    call printarray (rcollection)
    
    ! Remove the data
    call collct_deletevalue (rcollection,"MYIARRAY")
    
    ! Print statistics of the collection
    call collct_printStatistics (rcollection)

    ! Release the collection
    call collct_done (rcollection)
    
  end subroutine

end module
