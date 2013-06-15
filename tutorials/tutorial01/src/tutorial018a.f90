!##############################################################################
!# Tutorial 018a: Basic collection handling, quick-access arrays
!##############################################################################

module tutorial018a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use collection

  implicit none
  private
  
  public :: start_tutorial018a

contains

  ! ***************************************************************************

  subroutine printarray (rcollection)
  type(t_collection), intent(inout) :: rcollection
  
    integer :: i
    
    call output_line ("Maximum #integers: "//sys_siL(ubound(rcollection%IquickAccess,1),10))
    call output_line ("Maximum #doubles:  "//sys_siL(ubound(rcollection%DquickAccess,1),10))
    call output_line ("Maximum #strings:  "//sys_siL(ubound(rcollection%SquickAccess,1),10))
    
    ! Print the content
    do i=1,5
      call output_line (sys_si(rcollection%IquickAccess(i),5))
    end do

    do i=1,5
      call output_line (sys_sd(rcollection%DquickAccess(i),10))
    end do
  
    do i=1,3
      call output_line (trim(rcollection%SquickAccess(i)))
    end do

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial018a

    ! Declare some variables.
    type(t_collection) :: rcollection
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 018a")
    call output_separator (OU_SEP_MINUS)

    ! For using the quick-access arrays, the collection does not need to be
    ! initialised. Just use it.
    !
    ! Initialise the arrays.
    rcollection%IquickAccess(1:5) = (/ 5,4,3,2,1 /)
    
    rcollection%DquickAccess(1:5) = (/ 10.0_DP,9.0_DP,8.0_DP,7.0_DP,6.0_DP /)
    
    rcollection%SquickAccess(1) = "STRING1"
    rcollection%SquickAccess(2) = "STRING2"
    rcollection%SquickAccess(3) = "STRING3"
    
    ! Pass the collection to the subroutine
    call printarray (rcollection)
    
  end subroutine

end module
