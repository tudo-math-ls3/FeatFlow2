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
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 003d")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Clear and start the timer
    ! =================================

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)
    
    ! =================================
    ! Do some work
    ! =================================
    
    ! Wait two seconds
    do
      call stat_sampleTimer(rtimer,delapsedTime)
      if ( delapsedTime .ge. 2.0_DP) exit
    end do
    
    ! =================================
    ! Stop the timer
    ! =================================

    call stat_stopTimer (rtimer)
    
    call output_line ("Time for calculation in seconds: " // &
        trim(sys_sdL(dint(rtimer%delapsedReal),1)) )
    
  end subroutine

end module
