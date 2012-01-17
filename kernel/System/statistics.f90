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
!#  3.) stat_sampleTimer
!#      -> calculates the time that has elapsed since it was started without
!#         stopping the timer
!#
!#  4.) stat_stopTimer
!#      -> stops the timer object and calculates the time that has elapsed
!#         since it was started
!#
!#  5.) stat_sgetTime_byTimer
!#      -> Returns a string containing both CPU and real elapsed time
!#
!#  6.) stat_addTimers
!#      -> Computes the sum of the elapsed times the given timers hold
!#
!#  7.) stat_subTimers
!#      -> Computes the difference of the elapsed times (rtimer2 - rtimer1)
!#         the given timers hold.
!#
!#  8.) stat_rcloneTimer
!#      -> Clones the given timer
!#
!# Auxiliary routines:
!#
!#  1.) stat_calender_to_julian
!#      -> Reformat a standard date to a Julian date
!#
!#
!#  Usage of timers  \\
!# ----------------- \\
!# To calculate the computation time of another routine, one has to use
!# the t_timer structure. This can be done as follows:
!#
!# <code>
!#   ! Declare a timer structure
!#   TYPE(t_timer) :: rtimer
!#
!#   ! Initialise the timer by zero:
!#   CALL stat_clearTimer(rtimer)
!#
!#   ! Start the timer
!#   CALL stat_startTimer(rtimer)
!#
!#   ! Do some time-consuming work
!#   ...
!#
!#   ! Stop the timer
!#   CALL stat_stopTimer(rtimer)
!#
!#   ! Print the wall clock time that was necessary for the computation
!#   WRITE(*,*) "Time for computation: ",rtimer%delapsedReal
!# </code>
!#
!# Note that there may be some confusion in case, OpenMP is used! Use
!# rtimer%delapsedReal to get the real computation time of the program and
!# assure that you are the only one computing on the machine!
!#
!#  Short- and long-term timers
!# -----------------------------
!# The additional parameter CTYPE in stat_startTimer allows to configure
!# a timer to be a short-term or long-term timer. By default (if the
!# parameter is not specified), timers are both, short-term and long-term
!# timers. Nevertheless be aware, that using timers (especially standard
!# timers) cost processor cycles when using them. As a rule of thumb,
!#
!# a) when you are using smaller algorithms which need exact timings and
!#    do not want to waste time too much time in timer handling, use
!#    short-term timers by calling
!#
!# <code>
!#      CALL stat_startTimer(rtimer,STAT_TIMERSHORT)
!# </code>
!#
!#    These timers are much faster to handle than standard timers as they
!#    calculate processor cycles (at least in most cases, depending on the
!#    underlying processor architecture). Note that these timers are not
!#    designed for long-term time measurements as (depending on the underlying
!#    architecture) a timer may get an overflow!
!#
!# b) for larger algorithms which do not need exact timing, use either
!#    a long-term timer by calling
!#
!# <code>
!#      CALL stat_startTimer(rtimer,STAT_TIMERSHORT)
!# </code>
!#
!#    or a standard timer by calling
!#
!# <code>
!#      CALL stat_startTimer(rtimer,STAT_TIMERSTANDARD)
!# </code>
!#
!#    or
!#
!# <code>
!#      CALL stat_startTimer(rtimer).
!# </code>
!#
!#    A long-term timer calculates the time with a granularity of milliseconds.
!#    A standard timer automatically chooses the most accurate time measurement
!#    available -- i.e. for short time spans it uses a short-term timer and
!#    for longer time intervals a long-term timer. A pure long-term timer will
!#    calculate only the wall clock rtimer%delapsedReal, while a standard timer
!#    calculates both, wall clock (rtimer%delapsedReal) and CPU clock
!#    (rtimer%delapsedCPU).
!#
!#    Note that both, standard and long-term timers, are much more cost
!#    intensive in terms of processor cycles than short-term timers! They
!#    should be used to stop the time of larger algorithms or if the costs
!#    for an algorithm  cannot be guessed a-priori. A too intensive use of
!#    this type of timer may have an impact on the performance of the whole
!#    program!
!#
!# </purpose>
!##############################################################################

module statistics

!$use omp_lib
  use fsystem
  use genoutput

  implicit none
  
  private

!<constants>
  
  !<constantblock description="Timer types">
  
  ! Undefined timer, timer not running.
  integer, parameter, public :: STAT_TIMERNOTRUNNING = 0
  
  ! Short-term timer
  integer, parameter, public :: STAT_TIMERSHORT = 2**0
  
  ! Long-term timer
  integer, parameter, public :: STAT_TIMERLONG = 2**1
  
  ! Standard timer for both, short and long time time measurement
  integer, parameter, public :: STAT_TIMERSTANDARD = STAT_TIMERSHORT + STAT_TIMERLONG
  
  !</constantblock>

!</constants>

!<types>

!<typeblock>
  ! Timer object
  type t_timer
  
    ! Type of timer. A STAT_TIMERxxxx constant
    integer :: ctype = STAT_TIMERNOTRUNNING
    
    !<!-- Short time timer -->
  
    ! Short-term timer: Elapsed CPU time (clock cycles / frequency)
    ! warning: might be inaccurate for GPU code or
    ! parallel code
    real(DP) :: delapsedCPU = 0.0_DP
    
    ! Short-term timer: Elapsed real time (wall clock)
    real(DP) :: delapsedReal = 0.0_DP
    
    ! Short-term timer: Start CPU time (clock cycles / frequency)
    ! warning: might be inaccurate for GPU code or
    ! parallel code
    real(DP) :: dstartCPU = 0.0_DP

    ! Short-term timer: Start real time (wall clock)
    real(DP) :: dstartReal = 0.0_DP
    
    ! Value of sysclock counter during last call to stat_startTimer
    !  (to avoid floating point cancellation effects
    integer :: istartCount = 0
    
    !<!-- Long-term timer -->
    
    ! Long-term timer: Start date (wall clock); format is a Julian date.
    integer :: istartCPUdate = 0
    
    ! Long-term timer: Start time (wall clock) in ms relative to 0:00:00
    ! of the current day.
    integer :: istartCPUtime = 0
    
  end type t_timer
  
  public :: t_timer
!</typeblock>

!</types>

  public :: stat_clearTimer
  public :: stat_startTimer
  public :: stat_stopTimer
  public :: stat_sgetTime_byTimer
  public :: stat_addTimers
  public :: stat_subTimers
  public :: stat_rcloneTimer
  public :: stat_calender_to_julian
  public :: stat_sampleTimer
  
contains

! ***************************************************************************

!<subroutine>

  subroutine stat_init()

!<description>
  ! Initialises the statistics module.
!</description>

!</subroutine>

    ! Currently: nothing to do.

  end subroutine stat_init

! ***************************************************************************

!<subroutine>

  subroutine stat_done()

!<description>
  ! Cleans up the statistic module.
!</description>

!</subroutine>

    ! Currently: nothing to do.

  end subroutine stat_done

! ***************************************************************************

!<function>
  elemental function stat_calender_to_julian(year, month, day) result(ivalue)

!<description>
  ! Converts a standard date into a Julian date.
!</description>

!<input>
    ! Year
    integer, intent(in) :: year
    
    ! Month
    integer, intent(in) :: month
    
    ! Day
    integer, intent(in) :: day

    integer             :: ivalue
!</input>

!<result>
  ! Julian date.
!</result>

!</function>

    ivalue = day-32075+&
              1461*(year+4800+(month-14)/12)/4+&
              367*(month-2-((month-14)/12)*12)/12-&
              3*((year+4900+(month-14)/12)/100)/4
              
  end function stat_calender_to_julian

! ***************************************************************************
! Short-term timer routines
! ***************************************************************************

!<subroutine>
  subroutine stat_clearTimer(rtimer)
!<description>
  ! Clears the given timer object and sets its status to "not running".
!</description>
!<input>
    type (t_timer), intent(inout) :: rtimer
!</input>
!</subroutine>

    rtimer%delapsedCPU    = 0.0_DP
    rtimer%delapsedReal   = 0.0_DP
    rtimer%istartCount    = 0
    rtimer%dstartCPU      = 0.0_DP
    rtimer%dstartReal     = 0.0_DP
    rtimer%istartCPUdate  = 0.0_DP
    rtimer%istartCPUtime  = 0.0_DP
    rtimer%ctype = STAT_TIMERNOTRUNNING

  end subroutine stat_clearTimer

! ***************************************************************************

!<subroutine>
  subroutine stat_startTimer(rtimer,ctype)

!<description>
  ! Starts the given timer object. ctype specifies the type of the timer.
  ! Short-term timer are computationally more efficient than long-term
  ! timer (they are just faster when being used), but they do not allow
  ! to compute the total time for long-term time measurements on some
  ! computer systems.
!</description>

!<input>
    ! OPTIONAL: Type of the timer to start.
    ! STAT_TIMERSHORT: Start a short-term timer.
    ! STAT_TIMERLONG: Start a long-term timer.
    ! STAT_TIMERSTANDARD: Start both, short- and long-term timer.
    ! If not present, STAT_TIMERSTANDARD is assumed.
    integer, intent(in), optional :: ctype
!</input>

!<inputoutput>
    ! The timer to start.
    type (t_timer), intent(inout) :: rtimer
!</inputoutput>

!</subroutine>

    integer  :: icount, irate, icmax
    real(SP) :: dtime
    integer, dimension(8) :: itimeLong

    ! Set the type
    rtimer%ctype = STAT_TIMERSTANDARD
    if (present(ctype)) rtimer%ctype = ctype
    
    ! Start the short-term time measurement?
    if (iand(rtimer%ctype,STAT_TIMERSHORT) .ne. 0) then

      call system_clock(icount, irate, icmax)
      rtimer%dstartReal  = real(icount, DP) / real(irate, DP)
      rtimer%istartCount = icount

      call cpu_time(dtime)
      rtimer%dstartCPU = real(dtime, DP)
      
    end if
    
    ! Start the long-term time-measurement?
    if (iand(rtimer%ctype,STAT_TIMERLONG) .ne. 0) then
      
      call date_and_time(values=itimeLong)
      
      ! Get the Julian date
      rtimer%istartCPUdate = &
          stat_calender_to_julian(itimeLong(1), itimeLong(2), itimeLong(3))
          
      ! Get the time since 0:00:00 in ms.
      rtimer%istartCPUtime = itimeLong(5) * 60 * 60 * 1000 + &
                             itimeLong(6) * 60 * 1000 + &
                             itimeLong(7) * 1000 + &
                             itimeLong(8)
    end if

  end subroutine stat_startTimer

! ***************************************************************************

!<subroutine>

  subroutine stat_stopTimer(rtimer)
  
!<description>
  ! Stops the given timer object and accumulates the time since it was started
  ! with the elapsed time this object already holds. In other words:
  ! A sequence of start/stopTimer() calls results in the total accumulated elapsed
  ! time of each timed block of code, not in that of the last one.
!</description>

!<inputoutput>
    ! The timer to stop.
    type (t_timer), intent(inout) :: rtimer
!</inputoutput>

!</subroutine>

    real(DP), parameter :: dsecPerDay = 24.0_DP*60.0_DP*60.0_DP

    integer :: icount, irate, icmax, icpudate,icputime
    real(DP) :: dcurrentCpu, dcurrentReal, dtmp, delapsed
    real(SP) :: dtime
    integer, dimension(8) :: itimeLong

    if (rtimer%ctype .eq. STAT_TIMERNOTRUNNING) then
      call output_line ('Timer not running!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'stat_stoptimer')
      return
    end if

    delapsed = 0.0_DP

    ! Is that a short-term timer?
    if (iand(rtimer%ctype,STAT_TIMERSHORT) .ne. 0) then

      ! Ask the system clock
      call system_clock(icount, irate, icmax)
      dCurrentReal = real(icount, DP) / real(irate, DP)
      dtmp = dcurrentReal - rtimer%dstartReal

      if (icount .lt. rtimer%istartCount) then
        ! Add the maximum measurable timespan to that negative number to
        ! get the correct timespan.
        dtmp = dtmp + real(icmax,DP)/real(irate,DP)
      endif
      
      ! Save the wall clock time. This may be modified below...
      delapsed = max(dtmp, 0.0_DP)
      
      ! Get the CPU time.
      call cpu_time(dtime)
      dcurrentCPU = real(dtime, DP)
      rtimer%delapsedCPU = rtimer%delapsedCPU + (dcurrentCPU - rtimer%dstartCPU)
    end if
    
    ! Is that a long-term timer?
    if (iand(rtimer%ctype,STAT_TIMERLONG) .ne. 0) then
    
      call date_and_time(values=itimeLong)
      
      ! Get the Julian date
      iCPUdate = &
          stat_calender_to_julian(itimeLong(1), itimeLong(2), itimeLong(3))
          
      ! Get the time since 0:00:00 in ms.
      iCPUtime = itimeLong(5) * 60 * 60 * 1000 + &
                             itimeLong(6) * 60 * 1000 + &
                             itimeLong(7) * 1000 + &
                             itimeLong(8)
                             
      ! Calculate the actual time since the last start_timer:
      dtmp = real(iCPUdate - rtimer%istartCPUdate,dp) * dsecPerDay + &
             real(iCPUtime - rtimer%istartCPUtime,dp) * 0.001_DP
             
      ! If both, short- and long-term timer, are activated, take the
      ! wall time that better suits the correct time.
      if (iand(rtimer%ctype,STAT_TIMERSHORT) .ne. 0) then
      
        ! If the difference is more that let us say 100 sec, take
        ! the approximated time. This should also be the case if the
        ! timer from system_clock gets more than one overflow...
        ! If the difference is less, take the exact time.
        if (abs(dtmp-delapsed) .gt. 100.0_DP) then
          delapsed = max(dtmp, 0.0_DP)
        end if
      
      else
        ! That is our wall clock time.
        delapsed = max(dtmp, 0.0_DP)
      end if
      
    end if
    
    ! Sum up the wall clock time
    rtimer%delapsedReal = rtimer%delapsedReal + delapsed
    
    ! Switch off the timer.
    rtimer%dstartReal = 0.0_DP
    rtimer%dstartCPU = 0.0_DP
    rtimer%istartCount = 0
    rtimer%istartCPUdate = 0
    rtimer%istartCPUtime = 0
    rtimer%ctype = STAT_TIMERNOTRUNNING
      
  end subroutine stat_stopTimer
  
! ***************************************************************************

!<subroutine>

  subroutine stat_sampleTimer(rtimer,delapsedTime)
  
!<description>
  ! This routine calculates the wallclock time from starting the timer with
  ! stat_startTimer until now but does not stop the timer.
!</description>

!<input>
    ! The timer to stop.
    type (t_timer), intent(inout) :: rtimer
!</input>

!<output>
    ! Elapsed time since last start of the timer.
    real(DP), intent(out) :: delapsedTime
!</output>

!</subroutine>

    real(DP), parameter :: dsecPerDay = 24.0_DP*60.0_DP*60.0_DP

    integer :: icount, irate, icmax, icpudate,icputime
    real(DP) :: dcurrentReal, dtmp
    integer, dimension(8) :: itimeLong

    if (rtimer%ctype .eq. STAT_TIMERNOTRUNNING) then
      call output_line ('Timer not running!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'stat_stoptimer')
      return
    end if

    delapsedTime = 0.0_DP

    ! Is that a short-term timer?
    if (iand(rtimer%ctype,STAT_TIMERSHORT) .ne. 0) then

      ! Ask the system clock
      call system_clock(icount, irate, icmax)
      dCurrentReal = real(icount, DP) / real(irate, DP)
      dtmp = dcurrentReal - rtimer%dstartReal

      if (icount .lt. rtimer%istartCount) then
        ! Add the maximum measurable timespan to that negative number to
        ! get the correct timespan.
        dtmp = dtmp + real(icmax,DP)/real(irate,DP)
      endif
      
      ! Save the wall clock time. This may be modified below...
      delapsedTime = max(dtmp, 0.0_DP)
    end if
    
    ! Is that a long-term timer?
    if (iand(rtimer%ctype,STAT_TIMERLONG) .ne. 0) then
    
      call date_and_time(values=itimeLong)
      
      ! Get the Julian date
      iCPUdate = &
          stat_calender_to_julian(itimeLong(1), itimeLong(2), itimeLong(3))
          
      ! Get the time since 0:00:00 in ms.
      iCPUtime = itimeLong(5) * 60 * 60 * 1000 + &
                             itimeLong(6) * 60 * 1000 + &
                             itimeLong(7) * 1000 + &
                             itimeLong(8)
                             
      ! Calculate the actual time since the last start_timer:
      dtmp = real(iCPUdate - rtimer%istartCPUdate,dp) * dsecPerDay + &
             real(iCPUtime - rtimer%istartCPUtime,dp) * 0.001_DP
             
      ! If both, short- and long-term timer, are activated, take the
      ! wall time that better suits the correct time.
      if (iand(rtimer%ctype,STAT_TIMERSHORT) .ne. 0) then
      
        ! If the difference is more that let us say 100 sec, take
        ! the approximated time. This should also be the case if the
        ! timer from system_clock gets more than one overflow...
        ! If the difference is less, take the exact time.
        if (abs(dtmp-delapsedTime) .gt. 100.0_DP) then
          delapsedTime = max(dtmp, 0.0_DP)
        end if
      
      else
        ! That is our wall clock time.
        delapsedTime = max(dtmp, 0.0_DP)
      end if
      
    end if
    
  end subroutine stat_sampleTimer
  
  ! ***************************************************************************

!<function>
  function stat_sgetTime_byTimer(rtimer) result (stime)

!<description>
  ! Returns a string containing both CPU and real elapsed time.
  ! If the timer is running, the output will be nonsense. Without
  ! ENABLE_PARAMETER_CHECK, this is not even reported to improve performance.
!</description>

!<input>
    type (t_timer), intent(in) :: rtimer
!</input>

!<result>
    character(len=14) :: stime
!</result>

!</function>

    character(len=6) saux1, saux2

    saux1 = sys_sdL(rtimer%delapsedCPU, 2)
    saux2 = sys_sdL(rtimer%delapsedReal, 2)
    stime = adjustr(saux1) // " /" // adjustr(saux2)

  end function stat_sgetTime_byTimer

! ***************************************************************************

!<subroutine>
  subroutine stat_addTimers(rtimer1,rtimer2)

!<description>
  ! Computes the sum of the elapsed times the given timers hold and stores it in rtimer2.
  ! Both timers must be stopped, otherwise the result is nonsense. Without
  ! ENABLE_PARAMETER_CHECK, this is not even reported to improve performance.
!</description>

!<input>
    type (t_timer), intent(in) :: rtimer1
!</input>

!<inputoutput>
    type (t_timer), intent(inout) :: rtimer2
!</inputoutput>

!</subroutine>

    rtimer2%delapsedCPU  = rtimer1%delapsedCPU  + rtimer2%delapsedCPU
    rtimer2%delapsedReal = rtimer1%delapsedReal + rtimer2%delapsedReal
    rtimer2%dstartCPU    = 0.0_DP
    rtimer2%dstartReal   = 0.0_DP
    rtimer2%istartCPUdate = 0
    rtimer2%istartCPUtime = 0
    rtimer2%ctype = STAT_TIMERNOTRUNNING

  end subroutine stat_addTimers

! ***************************************************************************

!<subroutine>
  subroutine stat_subTimers(rtimer1,rtimer2)

!<description>
  ! Computes the difference of the elapsed times (rtimer2 = rtimer2 - rtimer1) the given
  ! timers represent.  Both timers must be stopped, otherwise the result is nonsense.
  ! Without ENABLE_PARAMETER_CHECK, this is not even reported to improve performance.
!</description>

!<input>
    type (t_timer), intent(in) :: rtimer1
!</input>

!<inputoutput>
    type (t_timer), intent(inout) :: rtimer2
!</inputoutput>

!</subroutine>

    rtimer2%delapsedCPU  = rtimer2%delapsedCPU  - rtimer1%delapsedCPU
    rtimer2%delapsedReal = rtimer2%delapsedReal - rtimer1%delapsedReal
    rtimer2%dstartCPU    = 0.0_DP
    rtimer2%dstartReal   = 0.0_DP
    rtimer2%istartCPUdate = 0
    rtimer2%istartCPUtime = 0
    rtimer2%ctype = STAT_TIMERNOTRUNNING

  end subroutine stat_subTimers

! ***************************************************************************

!<function>
  function stat_rcloneTimer(rtimer) result (rtimerResult)

!<description>
  ! Clones the given timer, which must be stopped.
  ! If the input timer is not stopped, otherwise the result is nonsense. Without
  ! ENABLE_PARAMETER_CHECK, this is not even reported to improve performance.
!</description>

!<input>
    type (t_timer), intent(in) :: rtimer
!</input>

!<result>
    type (t_timer) :: rtimerResult
!</result>

!</function>

    rtimerResult%delapsedCPU  = rtimer%delapsedCPU
    rtimerResult%delapsedReal = rtimer%delapsedReal
    rtimerResult%dstartCPU    = rtimer%dstartCPU
    rtimerResult%dstartReal   = rtimer%dstartReal
    rtimerResult%istartCPUdate = rtimer%istartCPUdate
    rtimerResult%istartCPUtime = rtimer%istartCPUtime
    rtimerResult%ctype = rtimer%ctype

  end function stat_rcloneTimer

end module
