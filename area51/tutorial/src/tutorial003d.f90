!##############################################################################
!# Tutorial 003d: Using the timer to measure time
!##############################################################################

module tutorial003d

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use statistics

  implicit none
  private
  
  public :: start_tutorial003d

contains

  ! ***************************************************************************

  subroutine start_tutorial003d

    ! Declare some variables
    type(t_timer) :: rtimer
    real(DP) :: delapsedTime
    
    integer :: i,d

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("Hello world. This is FEAT-2. Tutorial 003d")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Clear and start the timer
    ! =================================

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)
    
    ! =================================
    ! Do some work
    ! =================================
    do i = 1,100000000
      d = sqrt ( real(i,DP) )
    end do
    
    ! =================================
    ! Stop the timer
    ! =================================

    call stat_stopTimer (rtimer)
    
    call output_line ("Time for calculation in seconds: " // &
        trim(sys_sdL(dnint(rtimer%delapsedReal),1)) )
    
  end subroutine

end module
