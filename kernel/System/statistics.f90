!##############################################################################
!# ****************************************************************************
!# <name> statistics </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains system routines to construct timer objects
!# that can be used to calculate the time spent for code passages, functions
!# or similar units. The time can be measured in cpu cycles or clock cycles.
!# In the future 'statistic objects' can be included to calculate information
!# like MFLOPS and so on...
!#
!# The following routines can be found here
!#
!#  1.) stat_clearTimer
!#      -> clears the internal values of the timer object
!#
!#  2.) stat_startTimer
!#      -> starts the timer object
!#
!#  3.) stat_stopTimer
!#      -> stops the timer object and calculates the time that has elapsed since it was started
!#
!#  4.) stat_sgetTime2
!#      -> Returns a string containing both CPU and real elapsed time
!#
!#  5.) stat_addTimer
!#      -> Computes the sum of the elapsed times the given timers hold
!#
!#  6.) stat_subTimer
!#      -> Computes the difference of the elapsed times (rtimer2 - rtimer1) the given timers hold.
!#
!#  7.) stat_cloneTimer
!#      -> Clones the given timer
!#
!#  8.) stat_clearOp
!#      -> This routine clears the nop statistics 
!#
!#  9.) stat_clearTime
!#      -> This routine clears the time statistics 
!#
!# 10.) stat_clear
!#      -> This routine clears the statistic in rstat
!# 
!#
!#  ... (documentation incomplete)
!# </purpose>
!##############################################################################



module statistics

  use fsystem
  use genoutput


  implicit none
  
  !<constants>
  
  !<constantblock description="constants for special timer values">
  
  REAL(DP), PARAMETER :: STAT_TIMER_NOT_RUNNING = -1.0d0
  
  !</constantblock>
  
  !</constants>


!<types>

!<typeblock>
  ! Statistics about runtime and number of operation
  type t_stat

    !number of operations
    integer(I64) :: nflops

    !CPU time (clock cycles / frequency)
    real(DP) :: dtimeCpu

    !Real time (wall-clock)
    real(DP) :: dtimeReal 

    !number of memory loads (from RAM to CPU)
    integer(I64) :: nloads

    !number of memory stores (from CPU to RAM)
    integer(I64) :: nstores

  end type t_stat
!</typeblock>


!<typeblock>
  ! Timer object
  type t_timer

    !CPU time (clock cycles / frequency)
    !warning: might be inaccurate for GPU code.
    real(DP) :: delapsedCpu, dstartCpu

    !Real time (wall clock)
    real(DP) :: delapsedReal, dstartReal

  end type t_timer

!</typeblock>

!</types>

  ! maximal time value (for system_clock)
 ! real(DP) :: sys_dsimtime

!<globals>
  type(t_stat) :: sys_comminittime
  type(t_stat) :: sys_commgtime
  type(t_stat) :: sys_commitime
  type(t_stat) :: sys_commetime
 
  type(t_stat) :: stat_setBorderGTime

  type(t_stat) :: stat_matrixGTime
  type(t_stat) :: stat_matrixCorrTime
 
  type(t_stat) :: stat_smoothTime
  type(t_stat) :: stat_toolTime
  type(t_stat) :: stat_coarseTime
  type(t_stat) :: stat_scarcTime
  type(t_stat) :: stat_completeTime


  type(t_stat) :: stat_gridGenTime
  type(t_stat) :: stat_transferTime
  
  real(DP) :: stat_sizeSmoothData
  real(DP) :: stat_sizeMessageSent
  real(DP) :: stat_sizeSmoothingStepsPerf
  real(DP) :: stat_sizeTransferData

  type(t_stat) :: stat_assMatTime
  type(t_stat) :: stat_assRHSTime
!</globals>
  
  contains




!<subroutine>
  subroutine stat_clearTimer(rtimer)
! <description>
! Clears the given timer object and sets its status to "not running".
! </description>
! <input>
    type (t_timer) :: rtimer
! </input>
!</subroutine>

    rtimer%delapsedCpu = 0.0_DP
    rtimer%delapsedReal = 0.0_DP
    rtimer%dstartCpu = STAT_TIMER_NOT_RUNNING
    rtimer%dstartCpu = STAT_TIMER_NOT_RUNNING

  end subroutine stat_clearTimer


!<subroutine>
  subroutine stat_startTimer(rtimer)
! <description>
! Starts the given timer object.
! </description>
! <input>
    type (t_timer) :: rtimer
! </input>
!</subroutine>

    integer :: icount, irate, icmax

    call system_clock(icount, irate, icmax)
    rtimer%dstartReal = real(icount, DP) / real(irate, DP)

    call cpu_time(rtimer%dstartCpu)

  end subroutine stat_startTimer


!<subroutine>
  subroutine stat_stopTimer(rtimer)
! <description>
! Stops the given timer object and accumulates the time since it was started 
! with the elapsed time this object already holds. In other words: 
! A sequence of start/stopTimer() calls results in the total accumulated elapsed 
! time of each timed block of code, not in that of the last one.
! </description>
! <input>
    type (t_timer) :: rtimer
! </input>
!</subroutine>

    integer :: icount, irate, icmax
    real(DP) :: dcurrentCpu, dcurrentReal

    if (rtimer%dstartReal .ne. STAT_TIMER_NOT_RUNNING) then

      call system_clock(icount, irate, icmax)
      dCurrentReal = real(icount, DP) / real(irate, DP)
      rtimer%delapsedReal = rtimer%delapsedReal + (dcurrentReal - rtimer%dstartReal)
      rtimer%dstartReal = STAT_TIMER_NOT_RUNNING
      
      call cpu_time(dcurrentCpu)
      rtimer%delapsedCpu = rtimer%delapsedCpu + (dcurrentCpu - rtimer%dstartCpu)
      rtimer%dstartCpu = STAT_TIMER_NOT_RUNNING
      
    else

      call output_line(OU_CLASS_ERROR, "stat_stoptimer", &
       "You must not try to stop a timer that is not running.") 

    end if

  end subroutine stat_stopTimer
  
!<function>
  !TODO rename to stat_sgetTime when reworking of statistics is finished.
  character(len=14) function stat_sgetTime2(rtimer) result (stime)
! <description>
! Returns a string containing both CPU and real elapsed time.
! </description>
! <input>
    type (t_timer) :: rtimer
! </input>
! <result>
!   character(len=14) :: stime
! </result>
!</function>

    character(len=6) saux1, saux2

    if (rtimer%dstartCpu .eq. STAT_TIMER_NOT_RUNNING) then
      saux1 = sys_sdL(rtimer%delapsedCpu, 2)
      saux2 = sys_sdL(rtimer%delapsedReal, 2)
      stime = adjustr(saux1) // " /" // adjustr(saux2)
    else
       stime = "XXXXXXXXXX"
    end if
  end function stat_sgetTime2



!<function>
  type (t_timer) function stat_addTimer(rtimer1,rtimer2) result (rtimerResult)
! <description>
! Computes the sum of the elapsed times the given timers hold.
! </description>
! <input>
    type (t_timer) :: rtimer1, rtimer2
! </input>
! <result>
!   type (t_timer) :: rtimerResult
! </result>
!</function>


    if (rtimer1%dstartCpu .eq. STAT_TIMER_NOT_RUNNING .and. &
        rtimer2%dstartCpu .eq. STAT_TIMER_NOT_RUNNING  ) then
      rtimerResult%delapsedCpu  = rtimer1%delapsedCpu  + rtimer2%delapsedCpu
      rtimerResult%delapsedReal = rtimer1%delapsedReal + rtimer2%delapsedReal
      rtimerResult%dstartCpu = STAT_TIMER_NOT_RUNNING
      rtimerResult%dstartReal = STAT_TIMER_NOT_RUNNING
    else
      call output_line(OU_CLASS_ERROR, "stat_addtimer", "You must not try to add timers that are running.") 

    end if
  end function stat_addTimer


!<function>
  type (t_timer) function stat_subTimer(rtimer1,rtimer2) result (rtimerResult)
! <description>
! Computes the difference of the elapsed times (rtimer2 - rtimer1) the given timers hold.
! </description>
! <input>
    type (t_timer) :: rtimer1, rtimer2
! </input>
! <result>
!   type (t_timer) :: rtimerResult
! </result>
!</function>


    if (rtimer1%dstartCpu .eq. STAT_TIMER_NOT_RUNNING .and. &
        rtimer2%dstartCpu .eq. STAT_TIMER_NOT_RUNNING  ) then
      rtimerResult%delapsedCpu  = rtimer2%delapsedCpu  - rtimer1%delapsedCpu
      rtimerResult%delapsedReal = rtimer2%delapsedReal - rtimer1%delapsedReal
      rtimerResult%dstartCpu = STAT_TIMER_NOT_RUNNING
      rtimerResult%dstartReal = STAT_TIMER_NOT_RUNNING
    else
      call output_line(OU_CLASS_ERROR, "stat_subtimer", "You must not try to subtract timers that are running.") 

    end if
  end function stat_subTimer


!<function>
  type (t_timer) function stat_cloneTimer(rtimer) result (rtimerResult)
! <description>
! Clones the given timer.
! </description>
! <input>
    type (t_timer) :: rtimer
! </input>
! <result>
!   type (t_timer) :: rtimerResult
! </result>
!</function>


    if (rtimer%dstartCpu .eq. STAT_TIMER_NOT_RUNNING) then
      rtimerResult%delapsedCpu  = rtimer%delapsedCpu
      rtimerResult%delapsedReal = rtimer%delapsedReal
      rtimerResult%dstartCpu = rtimer%dstartCpu
      rtimerResult%dstartReal = rtimer%dstartReal
    else
      call output_line(OU_CLASS_ERROR, "stat_subtimer", "You must not try to clone a timer that is running.") 

    end if
  end function stat_cloneTimer

!************************************************************************

!<subroutine>
  subroutine stat_clearOp(rstat)

! <description>
! This routine clears the nop statistics in rstat
! </description>

! <input>
    type (t_stat) :: rstat
! </input>

!</subroutine>

    rstat%nflops  = 0_I64
    rstat%nloads  = 0_I64
    rstat%nstores = 0_I64

  end subroutine stat_clearOp

!************************************************************************

!<subroutine>
  subroutine stat_clearTime(rstat)

! <description>
! This routine clears the time statistics in rstat
! </description>

! <input>
    type (t_stat) :: rstat
! </input>

!</subroutine>

    rstat%dtimeCpu  = 0.0_DP
    rstat%dtimeReal = 0.0_DP

  end subroutine stat_clearTime

!************************************************************************

!<subroutine>
  subroutine stat_clear(rstat)

! <description>
! This routine clears the statistic in rstat
! </description>

! <input>
    type (t_stat) :: rstat
! </input>

!</subroutine>

    call stat_clearOp(rstat)
    call stat_clearTime(rstat)

  end subroutine stat_clear

!************************************************************************


end module
