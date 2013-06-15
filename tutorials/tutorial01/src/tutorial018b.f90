!##############################################################################
!# Tutorial 018b: Basic collection handling, linking collections
!##############################################################################

module tutorial018b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use collection

  implicit none
  private
  
  public :: start_tutorial018b

contains

  ! ***************************************************************************

  subroutine printdarray (rcollection)
  type(t_collection), intent(inout) :: rcollection

    integer :: i

    ! Print the data of the collection
    do i=1,5
      call output_line (sys_sd(rcollection%DquickAccess(i),10))
    end do

  end subroutine

  ! ***************************************************************************

  subroutine printiarray (rcollection)
  type(t_collection), intent(inout) :: rcollection

    integer :: i

    ! Print the data in the sub-collection
    call printdarray (rcollection%p_rnextCollection)
    
    ! Print the data in the main collection
    do i=1,5
      call output_line (sys_si(rcollection%IquickAccess(i),5))
    end do

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial018b

    ! Declare some variables.
    type(t_collection) :: rcollection
    type(t_collection), target :: rsubCollection
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 018b")
    call output_separator (OU_SEP_MINUS)

    ! Initialise the arrays.
    rcollection%IquickAccess(1:5) = (/ 5,4,3,2,1 /)
    
    rsubCollection%DquickAccess(1:5) = (/ 10.0_DP, 9.0_DP, 8.0_DP, 7.0_DP, 6.0_DP /)    
    
    ! rsubCollection is a "sub"-collection of rcollection
    rcollection%p_rnextCollection => rsubCollection
    
    ! Pass the collection to the subroutine
    call printiarray (rcollection)
    
  end subroutine

end module
