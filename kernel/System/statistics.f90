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
!#  5.) stat_addTimers
!#      -> Computes the sum of the elapsed times the given timers hold
!#
!#  6.) stat_subTimers
!#      -> Computes the difference of the elapsed times (rtimer2 - rtimer1) the given timers hold.
!#
!#  7.) stat_rcloneTimer
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
!# To calculate the computation time of another routine, one has to use
!# the t_timer structure. This can be done as follows:
!#
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
!#
!# Note that there may be some confusion in case, OpenMP is used! Use 
!# rtimer%delapsedReal to get the real computation time of the program and
!# assure that you are the only one computing on the machine!
!# 
!# </purpose>
!##############################################################################

MODULE statistics

  USE fsystem
  USE genoutput

  IMPLICIT NONE
  
!<constants>
  
!<constantblock description="constants for special timer values">
  
  REAL(DP), PARAMETER :: STAT_TIMER_NOT_RUNNING = -1.0d0
  
  INTEGER, PARAMETER :: STAT_MAXPERFCOUNTERTYPES = 100
  
  INTEGER, PARAMETER :: STAT_MAXPERFCOUNTERS = 100
  
  INTEGER, PARAMETER :: STAT_PERFCOUNTER_SCARC_DYN = 100
  
  INTEGER, PARAMETER :: STAT_PERFCOUNTER_MULTIDIM_DYN = 200
  
  INTEGER, PARAMETER :: STAT_SCARCPERFCOUNTERS = 2
  
  INTEGER, PARAMETER :: STAT_MULTIDIMPERFCOUNTERS = 3
  
  INTEGER, PARAMETER :: STAT_STATICPERFCOUNTERS = 1
  
  INTEGER, PARAMETER :: STAT_SETINSITU = 0
  
  INTEGER, PARAMETER :: STAT_COPY =       1
  INTEGER, PARAMETER :: STAT_SCALEDCOPY = 2
  INTEGER, PARAMETER :: STAT_XPY =        3
  INTEGER, PARAMETER :: STAT_NXPY =       4
  INTEGER, PARAMETER :: STAT_AXPY =       5
  INTEGER, PARAMETER :: STAT_AXPYV =      6
  INTEGER, PARAMETER :: STAT_AXPBY =      7
  INTEGER, PARAMETER :: STAT_DOT =        8
  INTEGER, PARAMETER :: STAT_L1NORM =     9
  INTEGER, PARAMETER :: STAT_L2NORM =     10
  INTEGER, PARAMETER :: STAT_L2NORMSQUARED = 11
  INTEGER, PARAMETER :: STAT_MV_BAND_Q1_2D = 12
  
  INTEGER, PARAMETER :: STAT_MVONLY_BAND_Q1_2D=13
  INTEGER, PARAMETER :: STAT_MV_BAND_GENERIC_SBAND_2D=14
  INTEGER, PARAMETER :: STAT_PRECJAC_Q1_2D=15
  INTEGER, PARAMETER :: STAT_PRECGS_BAND_Q1_2D=16
  INTEGER, PARAMETER :: STAT_PRECTRIDI_BAND_Q1_2D=17
  INTEGER, PARAMETER :: STAT_PRECTRIGS_BAND_Q1_2D=18
  INTEGER, PARAMETER :: STAT_PRECILU_BAND_Q1_2D=19
  INTEGER, PARAMETER :: STAT_UMFPACK_2D=20
  INTEGER, PARAMETER :: STAT_PROLONGATE_Q1_2D=21
  INTEGER, PARAMETER :: STAT_RESTRICT_Q1_2D=22
  INTEGER, PARAMETER :: STAT_PROLONGATE_Q2_2D=23
  INTEGER, PARAMETER :: STAT_RESTRICT_Q2_2D=24
  INTEGER, PARAMETER :: STAT_PROLONGATE_Q2L_2D=25
  INTEGER, PARAMETER :: STAT_RESTRICT_Q2L_2D=26    
  INTEGER, PARAMETER :: STAT_NONTRIVIALFLOP=16
  INTEGER, PARAMETER :: STAT_PERFCOUNTER_SCARC_COMPLETE=1
  INTEGER, PARAMETER :: STAT_PERFCOUNTER_MDIM_COMPLETE=2

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


!<type>
  ! Timer object
  TYPE t_timer
    ! elapsed CPU time (clock cycles / frequency)
    ! warning: might be inaccurate for GPU code or
    ! parallel code
    real(DP) :: delapsedCPU != 0.0_DP
    ! start CPU time (clock cycles / frequency)
    ! warning: might be inaccurate for GPU code or
    ! parallel code
    real(DP) :: dstartCPU != STAT_TIMER_NOT_RUNNING

    ! elapsed real time (wall clock)
    real(DP) :: delapsedReal != 0.0_DP
    ! start real time (wall clock)
    real(DP) :: dstartReal != STAT_TIMER_NOT_RUNNING
    ! value of sysclock counter during last call to stat_startTimer
    !  (to avoid floating point cancellation effects
    integer :: istartCount
  END TYPE t_timer
!</type>

!<type>
  ! Performance counter object (combined timings, flops and memory transfers)
  ! This implementation is recursion-safe in the sense that it maintains
  ! a recursion counter and can be enabled and disabled as often as
  ! the code requires. It is finally deactivated once its recursion counter is
  ! back to zero.
  type t_perfCounter

    !timer
    type(t_timer) :: rtimer

    !number of floating point operations
    integer(I64)  :: nflops

    !number of memory loads (from RAM to CPU)
    integer(I64)  :: nloads

    !number of memory stores (from CPU to RAM)
    integer(I64)  :: nstores

    ! recursion counter (only used for static perfcounters)
    ! TODO maybe merge with inext below, which is only used for dynamic ones)
    integer       :: irecursionGuard

  end type t_perfCounter
!</type>
  
!<type>
  ! Internal type to maintain all arrays for performance counters easily.
  ! Users should typically not use this, instead, use the getter and 
  ! pretty-printer functions below.
  type t_perfCounterData

    ! Array for automatically collected statistics,
    ! indexed by globally defined PERFCOUNTER_BLABLA IDs.
    type(t_perfcounter), dimension(STAT_MAXPERFCOUNTERS) :: Rpc

    ! Bitfield to keep track of enabled and disabled perfcounters,
    ! indexed by the same globally defined PERFCOUNTER_BLABLA IDs.
    integer(I64) :: Ibits

    ! Pointer to next available perfcounter in case of dynamically maintained
    ! STAT_SCARCPERFCOUNTERS or STAT_MULTIDIMPERFCOUNTERS. Unused otherwise.
    integer :: inext

  end type t_perfCounterData
  private :: t_perfCounterData
!</type>
!</types>

  
!<privatevars>
  ! performance counter arrays
  type(t_perfCounterData), dimension(STAT_MAXPERFCOUNTERTYPES), private, target :: stat_RperfCounters
  ! The Fortran system_clock timer, like all integer timers, has a cycle time of 
  ! real(max)/real(rate) seconds. After max clock cycles the clock will start 
  ! counting again from zero so if you start getting negative times for sections 
  ! of your code you have probably exceeded the maximum amount of time that the 
  ! call can handle, and you have to add this value. All t_timer objects perform
  ! this automatically.
  real(DP), private :: stat_dtimeMax
!</privatevars>  

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

  ! new statistics: matrix assembly time
  type(t_timer) :: stat_matrixAssemblyTime
  ! new statistics: RHS assembly time
  type(t_timer) :: stat_rhsAssemblyTime
  ! new statistics: total time for solving coarse grid problems
  type(t_timer) :: stat_coarseTime2
  
  !<interfaces>
  ! description1
  interface stat_clearPerfCounter
    module procedure stat_clearPerfCounter_byID
    module procedure stat_clearPerfCounter_byRef
  end interface

  ! description2
  interface stat_rclonePerfCounter
    module procedure stat_rclonePerfCounter_byID
    module procedure stat_rclonePerfCounter_byRef
  end interface

  ! description3
  interface stat_sgetMFlopSec
    module procedure stat_sgetMFlopSec_byID
    module procedure stat_sgetMFlopSec_byRef
  end interface

  ! description4
  interface stat_sgetTime
    module procedure stat_sgetTime_byPCRef
    module procedure stat_sgetTime_byPCID
    module procedure stat_sgetTime_byTimer
  end interface

  ! description5
  interface stat_sgetTimeRatio
    module procedure stat_sgetTimeRatio_byPCRef
  end interface

  ! description6
  interface stat_ngetFlops
    module procedure stat_ngetFlops_byRef
    module procedure stat_ngetFlops_byID
  end interface
!</interfaces>

! corresponding private subroutines
  private :: stat_clearPerfCounter_byID
  private :: stat_clearPerfCounter_byRef
  private :: stat_rclonePerfCounter_byID
  private :: stat_rclonePerfCounter_byRef
  private :: stat_sgetMFlopSec_byID
  private :: stat_sgetMFlopSec_byRef
  private :: stat_sgetTime_byPCRef
  private :: stat_sgetTime_byPCID
  private :: stat_sgetTime_byTimer
  private :: stat_sgetTimeRatio_byPCRef
  private :: stat_ngetFlops_byRef
  private :: stat_ngetFlops_byID
  private :: stat_getCost

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

    rtimer%delapsedCPU    = 0.0_DP
    rtimer%delapsedReal   = 0.0_DP
    rtimer%istartCount    = 0
    rtimer%dstartCPU      = STAT_TIMER_NOT_RUNNING
    rtimer%dstartReal     = STAT_TIMER_NOT_RUNNING

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

    integer  :: icount, irate, icmax
    real(SP) :: dtime

    call system_clock(icount, irate, icmax)
    rtimer%dstartReal  = real(icount, DP) / real(irate, DP)
    rtimer%istartCount = icount

    call cpu_time(dtime)
    rtimer%dstartCPU = real(dtime, DP)

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
    real(DP) :: dcurrentCpu, dcurrentReal, dtmp
    real(SP) :: dtime

    if (rtimer%dstartReal .ne. STAT_TIMER_NOT_RUNNING) then

    call system_clock(icount, irate, icmax)
    dCurrentReal = real(icount, DP) / real(irate, DP)
    dtmp = dcurrentReal - rtimer%dstartReal

    if (icount .lt. rtimer%istartCount) then
      dtmp = dtmp + stat_dtimeMax
    endif
    rtimer%delapsedReal = rtimer%delapsedReal + max(dtmp, 0.0_DP)
    rtimer%dstartReal = STAT_TIMER_NOT_RUNNING
    rtimer%istartCount = 0

    call cpu_time(dtime)
    dcurrentCPU = real(dtime, DP)
    rtimer%delapsedCPU = rtimer%delapsedCPU + (dcurrentCPU - rtimer%dstartCPU)
    rtimer%dstartCPU = STAT_TIMER_NOT_RUNNING
      
    else

      call output_line(OU_CLASS_ERROR, "stat_stoptimer", &
       "You must not try to stop a timer that is not running.") 

    end if

  end subroutine stat_stopTimer
  
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
  !<errors>
  ! ERR_STAT_TIMERNOTRUNNING if timer is running (warning)
  !</errors>
!</function>

    character(len=6) saux1, saux2



    saux1 = sys_sdL(rtimer%delapsedCPU, 2)
    saux2 = sys_sdL(rtimer%delapsedReal, 2)
    stime = adjustr(saux1) // " /" // adjustr(saux2)


  end function stat_sgetTime_byTimer



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
  !<inoutput>
    type (t_timer) :: rtimer2
  !</inoutput>
  !<errors>
  ! ERR_STAT_TIMERRUNNING if any of the two timers is running (warning)
  !</errors>
!</subroutine>



    rtimer2%delapsedCPU  = rtimer1%delapsedCPU  + rtimer2%delapsedCPU
    rtimer2%delapsedReal = rtimer1%delapsedReal + rtimer2%delapsedReal
    rtimer2%dstartCPU    = STAT_TIMER_NOT_RUNNING
    rtimer2%dstartReal   = STAT_TIMER_NOT_RUNNING


  end subroutine stat_addTimers



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
  !<inoutput>
    type (t_timer) :: rtimer2
  !</inoutput>
  !<errors>
  ! ERR_STAT_TIMERRUNNING if any of the two timers is running (warning)
  !</errors>
!</subroutine>



    rtimer2%delapsedCPU  = rtimer2%delapsedCPU  - rtimer1%delapsedCPU
    rtimer2%delapsedReal = rtimer2%delapsedReal - rtimer1%delapsedReal
    rtimer2%dstartCPU    = STAT_TIMER_NOT_RUNNING
    rtimer2%dstartReal   = STAT_TIMER_NOT_RUNNING


  end subroutine stat_subTimers

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
  !<errors>
  ! ERR_STAT_TIMERRUNNING if the timer is running (warning)
  !</errors>
!</function>


    rtimerResult%delapsedCPU  = rtimer%delapsedCPU
    rtimerResult%delapsedReal = rtimer%delapsedReal
    rtimerResult%dstartCPU    = rtimer%dstartCPU
    rtimerResult%dstartReal   = rtimer%dstartReal


  end function stat_rcloneTimer

!************************************************************************

! ****************************************************************************************
!        !!!!!!!!!!!!!!!!!!!!!!!!! PERFCOUNTER STUFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ****************************************************************************************


!<subroutine>
  subroutine stat_init()
  !<description>
  ! Clears all performance counters unconditionally. Side effects are expected if
  ! this routine is called from anywhere except msrv_feastinit().
  !</description>
  !<errors>
  ! ERR_STAT_SYSTEMCLOCK if system_clock() seems unsupported (warning)
  !</errors>
!</subroutine>

    ! local loop variables
    integer :: icount, irate, icmax

    ! compute maximam measurable timespan
    call system_clock(icount, irate, icmax)

    ! compute maximam measurable timespan and do some statistics info output
    stat_dtimeMax = real(icmax,DP)/real(irate,DP)
!    call output_line(OL_MSG, "stat_init", &
!                     "Longest measurable timespam : " // &
!                     trim(sys_sdEL(stat_dtimeMax,2)) // " s.")
!    call output_line(OL_MSG, "stat_init", &
!                     "Shortest measurable timespam: " // &
!                     trim(sys_sdEL(1.0_DP/real(irate,DP),2)) // " s.")


  end subroutine stat_init


! ****************************************************************************************


!<subroutine>
  subroutine stat_clearPerfCounter_byID(cperfCounterID)
  !<description>
  ! Clears the given performance counter. 
  ! If the ID represents a static perfcounter, it is cleared only if it is not currently 
  ! counting any operations. Otherwise, an error message is displayed.
  ! If the given ID represents a dynamic performance counter 
  ! (STAT_PERFCOUNTER_SCARC_DYN  or STAT_PERFCOUNTER_MULTIDIM_DYN), 
  ! the innermost counter (counting from calls to stat_enable()) is cleared if it is not 
  ! currently counting any operations. Otherwise, an error message is displayed.
  !</description>
  !<input>
    ! one of the STAT_PERFCOUNTER_BLABLA IDs
    integer, intent(in) :: cperfCounterID
  !</input>
  !<errors>
  ! ERR_STAT_INVALID_STATICCOUNTER if perfcounter is unknown (critical)
  ! ERR_STAT_PERFCOUNTER_ACTIVE if perfcounter is active (critical)
  ! ERR_STAT_INVALID_DYNCOUNTER if perfcounter is unknown (critical)
  !</errors>
  !</subroutine>

    ! local vars
    integer :: iarray,icurrent

    ! figure out which array to work on
    if (cperfCounterID .eq. STAT_PERFCOUNTER_SCARC_DYN) then
      iarray = STAT_SCARCPERFCOUNTERS
    else if (cperfCounterID .eq. STAT_PERFCOUNTER_MULTIDIM_DYN) then
      iarray = STAT_MULTIDIMPERFCOUNTERS
    else 
      iarray = STAT_STATICPERFCOUNTERS  
    endif

    ! finally, clear the perfcounter
    if (iarray .eq. STAT_STATICPERFCOUNTERS) then
      icurrent = cperfCounterID
    else     
      icurrent=stat_RperfCounters(iarray)%inext-1
    endif

    ! error checking not necessary, will be done in stat_clearTimer()
    call stat_clearTimer(stat_RperfCounters(iarray)%Rpc(icurrent)%rtimer)
    stat_RperfCounters(iarray)%Rpc(icurrent)%nflops  = 0
    stat_RperfCounters(iarray)%Rpc(icurrent)%nloads  = 0
    stat_RperfCounters(iarray)%Rpc(icurrent)%nstores = 0

  end subroutine stat_clearPerfCounter_byID


! ****************************************************************************************


!<subroutine>
  subroutine stat_clearPerfCounter_byRef(rperfCounter)
  !<description>
  ! Clears the given performance counter. 
  !</description>
  !<inoutput>
    ! perfcounter to clear
    type(t_perfCounter) :: rperfCounter
  !</inoutput>
  !<errors>
  !  none
  !</errors>
  !</subroutine>

    ! error checking not necessary, will be done in stat_clearTimer()
    call stat_clearTimer(rperfCounter%rtimer)
    rperfCounter%nflops  = 0
    rperfCounter%nloads  = 0
    rperfCounter%nstores = 0

  end subroutine stat_clearPerfCounter_byRef


! ****************************************************************************************


!<function>
  function stat_rclonePerfCounter_byRef(rpc) result (rclone)
  !<description>
  ! Clones the given perfcounter.
  !</description>
  !<input>
    type(t_perfCounter), intent(in) :: rpc
  !</input>
  !<result>
    type(t_perfCounter) :: rclone
  !</result>
!</function>


    ! error checking not necessary, will be done in stat_cloneTimer()
    rclone%nflops          = rpc%nflops
    rclone%nloads          = rpc%nloads
    rclone%nstores         = rpc%nstores
    rclone%rtimer          = stat_rcloneTimer(rpc%rtimer)
    rclone%irecursionGuard = rpc%irecursionGuard


  end function stat_rclonePerfCounter_byRef


! ****************************************************************************************


!<function>
  function stat_rclonePerfCounter_byID(cperfCounterID) result (rclone)
  !<description>
  ! Clones the given static perfcounter. If the ID represents a dynamic perfcounter,
  ! a critical error is generated.
  !</description>
  !<input>
    integer, intent(in) :: cperfCounterID
  !</input>
  !<result>
    type(t_perfCounter) :: rclone
  !</result>
  !<errors>
  ! ERR_STAT_ILLEGAL_DYNCOUNTER if this routine is called with a dynamic counter
  ! ERR_STAT_INVALID_STATICCOUNTER if the given static perfcounter does not exist
  !</errors>
!</function>

    ! error checking not necessary, will be done in stat_cloneTimer()
    rclone%nflops          = stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                                 Rpc(cperfCounterID)%nflops
    rclone%nloads          = stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                                 Rpc(cperfCounterID)%nloads
    rclone%nstores         = stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                                 Rpc(cperfCounterID)%nstores
    rclone%rtimer          = stat_rcloneTimer(stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                                 Rpc(cperfCounterID)%rtimer)
    rclone%irecursionGuard = stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                                 Rpc(cperfCounterID)%irecursionGuard

  end function stat_rclonePerfCounter_byID


! ****************************************************************************************


!<subroutine>
  subroutine stat_addPerfCounters(rpc1, rpc2) 
  !<description>
  ! Adds the given perfcounters: rpc2 = rpc1 + rpc2
  !</description>
  !<input>
    type(t_perfCounter), intent(in) :: rpc1
  !</input>
  !<inoutput>
    type(t_perfCounter) :: rpc2
  !</inoutput>
  !<errors>
  ! none
  !</errors>
!</subroutine>

    ! error checking not necessary, will be done in stat_addTimers()
    call stat_addTimers(rpc1%rtimer, rpc2%rtimer)
    rpc2%nflops  = rpc1%nflops  + rpc2%nflops
    rpc2%nloads  = rpc1%nloads  + rpc2%nloads
    rpc2%nstores = rpc1%nstores + rpc2%nstores

  end subroutine stat_addPerfCounters


! ****************************************************************************************


!<subroutine>
  subroutine stat_enable(cperfCounterID)
  !<description>
  ! Enables the performance counter identified by cperfCounterID.
  !</description>
  !<input>
    ! one of the STAT_PERFCOUNTER_BLABLA IDs
    integer, intent(in) :: cperfCounterID
  !</input>
  !<errors>
  ! ERR_STAT_INVALID_STATICCOUNTER if perfcounter is unknown (critical)
  ! ERR_STAT_OVERFLOW if no dynamic perfcounters are left (critical)
  !</errors>
!</subroutine>

    ! local vars
    integer :: iarray,icurrent


    ! figure out which array to work on
    if (cperfCounterID .eq. STAT_PERFCOUNTER_SCARC_DYN) then
      iarray = STAT_SCARCPERFCOUNTERS
    else if (cperfCounterID .eq. STAT_PERFCOUNTER_MULTIDIM_DYN) then
      iarray = STAT_MULTIDIMPERFCOUNTERS
    else 
      iarray = STAT_STATICPERFCOUNTERS  
    endif

    ! finally, enable the counter
    if (iarray .eq. STAT_STATICPERFCOUNTERS) then

      ! check if this counter has to be enabled
      if (stat_RperfCounters(iarray)%Rpc(cperfCounterID)%irecursionGuard .eq. 0) then
        stat_RperfCounters(iarray)%Ibits = stat_RperfCounters(iarray)%Ibits + 2**cperfCounterID
        call stat_startTimer(stat_RperfCounters(iarray)%Rpc(cperfCounterID)%rtimer)
      endif
      ! update recursion counter
      stat_RperfCounters(iarray)%Rpc(cperfCounterID)%irecursionGuard = &
        stat_RperfCounters(iarray)%Rpc(cperfCounterID)%irecursionGuard + 1

    else

      icurrent = stat_RperfCounters(iarray)%inext
      call stat_startTimer(stat_RperfCounters(iarray)%Rpc(icurrent)%rtimer)
      stat_RperfCounters(iarray)%inext = stat_RperfCounters(iarray)%inext + 1

    endif

  end subroutine stat_enable


! ****************************************************************************************


!<subroutine>
  subroutine stat_disable(cperfCounterID, rpc)
  !<description>
  ! Disables the performance counter identified by cperfCounterID.
  ! For dynamic perf counters, the information is lost (and cleared) after disabling it. 
  ! If the optional argument is present, a clone of the counter is returned.
  !</description>
  !<input>
    ! one of the STAT_PERFCOUNTER_BLABLA IDs
    integer, intent(in) :: cperfCounterID
  !</input>
  !<output>
    ! a clone of the perfcounter that is disabled
    type(t_perfCounter), optional :: rpc
  !</output>
  !<errors>
  ! ERR_STAT_INVALID_STATICCOUNTER if perfcounter is unknown (critical)
  ! ERR_STAT_INVALID_DYNCOUNTER if perfcounter is unknown (critical)
  !</errors>
!</subroutine>

    ! local vars
    integer :: iarray,icurrent


    ! figure out which array to work on
    if (cperfCounterID .eq. STAT_PERFCOUNTER_SCARC_DYN) then
      iarray = STAT_SCARCPERFCOUNTERS
    else if (cperfCounterID .eq. STAT_PERFCOUNTER_MULTIDIM_DYN) then
      iarray = STAT_MULTIDIMPERFCOUNTERS
    else 
      iarray = STAT_STATICPERFCOUNTERS  
    endif

    ! finally, disable the counter
    if (iarray .eq. STAT_STATICPERFCOUNTERS) then

      ! update recursion counter
      stat_RperfCounters(iarray)%Rpc(cperfCounterID)%irecursionGuard = &
        stat_RperfCounters(iarray)%Rpc(cperfCounterID)%irecursionGuard - 1
      ! check if this counter can be disabled
      if (stat_RperfCounters(iarray)%Rpc(cperfCounterID)%irecursionGuard .eq. 0) then
        stat_RperfCounters(iarray)%Ibits = stat_RperfCounters(iarray)%Ibits - 2**cperfCounterID
        call stat_stopTimer(stat_RperfCounters(iarray)%Rpc(cperfCounterID)%rtimer)
      endif
      ! clone if necessary
      if (present(rpc)) then
        rpc = stat_rclonePerfCounter(stat_RperfCounters(iarray)%Rpc(cperfCounterID))
      endif

    else

      icurrent = stat_RperfCounters(iarray)%inext-1
      ! stop timer
      call stat_stopTimer(stat_RperfCounters(iarray)%Rpc(icurrent)%rtimer)
      ! clone if necessary
      if (present(rpc)) then
        rpc = stat_rclonePerfCounter(stat_RperfCounters(iarray)%Rpc(icurrent))
      endif
      ! remove data
      call stat_clearTimer(stat_RperfCounters(iarray)%Rpc(icurrent)%rtimer)
      stat_RperfCounters(iarray)%Rpc(icurrent)%nflops  = 0
      stat_RperfCounters(iarray)%Rpc(icurrent)%nloads  = 0
      stat_RperfCounters(iarray)%Rpc(icurrent)%nstores = 0
      ! and remove counter
      stat_RperfCounters(iarray)%inext = stat_RperfCounters(iarray)%inext - 1
    endif

  end subroutine stat_disable


! ****************************************************************************************


!<subroutine>
  subroutine stat_getCost(ctype, n, cprec, bconst, nflops, nloads, nstores)
  !<description>
  ! This routine returns for a specific operation and number of unknowns
  ! the number of floating point operations, and the amount of 
  ! floating point data moved to (nloads) and from (nstores) the CPU.
  ! Our usual abstraction model applies.
  ! </description>
  !<input>
    ! identifier / opcode
    integer, intent(in) :: ctype
    ! vector length
    integer, intent(in) :: n
    ! precision (STAT_SIZEOFFLOAT or STAT_SIZEOFDOUBLE)
    integer, intent(in) :: cprec
    ! constmatmul or varmatmul
    logical, intent(in) :: bconst
  !</input>
  !<output>
    ! flops for this operation
    integer, intent(out) :: nflops
    ! loads for this operation
    integer, intent(out) :: nloads
    ! stores for this operation
    integer, intent(out) :: nstores
  !</output>
  !<errors>
  ! ERR_STAT_INVALID_OPCODE if ctype is unknown (critical)
  !</errors>
!</subroutine>

    select case (ctype)

    !description syntax: 
    !    greek: scalars
    !    roman: vectors of length n
    !    capital roman: matrices (band length n, otherwise hard-coded constant)

    case(STAT_SETINSITU) 
      !x = alpha      (corresponding HPLA: dset)
      nflops  = 0 
      nloads  = cprec*1
      nstores = cprec*n

    case(STAT_COPY) 
      !y = x          (corresponding BLAS: dcopy)
      nflops  = 0
      nloads  = cprec*n
      nstores = cprec*n

    case(STAT_SCALEDCOPY) 
      !y = alpha*x  (corresponding HPLA: DCOPYS, corresponding BLAS: DSCAL [in situ])
      nflops  = n 
      nloads  = cprec*(n+1)
      nstores = cprec*n

    case(STAT_XPY)  
      !b = a+b, a-b, b-a, a*b  (corresponding HPLA: DXPY, DXSY, DXSY2,DXMY)
      nflops  = n
      nloads  = cprec*(2*n)
      nstores = cprec*n

    case(STAT_AXPYV)  
      !b = a+b*c
      nflops  = 2*n
      nloads  = cprec*(3*n)
      nstores = cprec*n

    case(STAT_NXPY)  
      !b = -(a+b)  (corresponding HPLA: DNXPY)
      nflops  = n+1
      nloads  = cprec*(2*n + 1)
      nstores = cprec*n

    case(STAT_AXPY)  
      !y = y + alpha*x  (corresponding BLAS: DAXPY)
      nflops  = 2*n
      nloads  = cprec*(2*n + 1)
      nstores = cprec*n

    case(STAT_AXPBY) 
      !b = alpha*a + beta*b (corresponding HPLA: DAXPBY)
      nflops  = 3*n
      nloads  = cprec*(2*n + 2)
      nstores = cprec*n

    case(STAT_DOT)  
      !gamma = a dot b (corresponding BLAS: ddot)
      nflops  = 2*n
      nloads  = cprec*(2*n)
      nstores = cprec*1

    case(STAT_L2NORM) 
      !gamma = sqrt(a dot a) (corresponding BLAS: DNRM2)
      nflops  = 2*n + STAT_NONTRIVIALFLOP
      nloads  = cprec*n 
      nstores = cprec*1

    case(STAT_L2NORMSQUARED) 
      !gamma = a dot a   (corresponding HPLA: DNRM22)
      nflops  = 2*n
      nloads  = cprec*n
      nstores = cprec*1

    case(STAT_L1NORM) 
      !gamma = sum_i x(i)  abs is for free
      nflops  = n
      nloads  = cprec*n
      nstores = cprec*1

    case(STAT_MV_BAND_Q1_2D)
      ! banded and const banded matrix-vector multiplication
      ! y = alpha*A*x + beta*y
      ! corresponding routine somewhere in matrix.f90
      if (bconst) then
        ! All additional ops and loads for special const treatment 
        ! are not counted, only the load for the matrix is reduced 
        ! from 9n to 9 according to our perf model.
        nflops  = 20*n  
        nloads  = cprec*(10*n)
        nstores = cprec*n 
      else
        nflops  = 20*n  
        nloads  = cprec*(19*n) 
        nstores = cprec*n
      endif

    case(STAT_MVONLY_BAND_Q1_2D)
      ! banded and const banded matrix-vector multiplication
      ! y = A*x
      ! corresponding routine somewhere in matrix.f90
      if (bconst) then
        ! All additional ops and loads for special const treatment 
        ! are not counted, only the load for the matrix is reduced 
        ! from 9n to 9 according to our perf model.
        nflops  = 17*n  
        nloads  = cprec*(9*n)
        nstores = cprec*n 
      else
        nflops  = 17*n  
        nloads  = cprec*(18*n) 
        nstores = cprec*n
      endif

    case(STAT_MV_BAND_GENERIC_SBAND_2D)
      ! multiplication of a single matrix band with a vector,
      ! y=y+alpha*band*x 
      ! corresponding routine somewhere in matrix.f90
      if (bconst) then
        nflops  = 3*n  
        nloads  = cprec*(2*n + 1)
        nstores = cprec*n 
     else
        nflops  = 3*n  
        nloads  = cprec*(3*n)
        nstores = cprec*n 
     endif

    case(STAT_PRECJAC_Q1_2D) 
      ! x = x + omega*c*d (corresponding HPLA: DJAC for cdefcorr.eq.0)
      ! Note: DJAC for cdefcorr.ne.0 is STAT_XPY
      if (bconst) then
        nflops  = 2*n + 1
        nloads  = cprec*(2*n + 2)
        nstores = cprec*n
      else
        nflops  = 3*n     
        nloads  = cprec*(3*n + 1)
        nstores = cprec*n
      endif

    case(STAT_PRECGS_BAND_Q1_2D)
      ! one Gauss-Seidel preconditioning step (corresponding routine: ib217.f)
      nflops  = 9*2*n + n
      nloads  = cprec*( (9*3*n) + (2*n) )
      nstores = cprec*n

    !MIGHT_BE_BUG: statistics is missing
    case(STAT_PRECTRIDI_BAND_Q1_2D) 
      ! prec_perform_tridi -> PRECS_LINE_E11_2D in sbblas
      nflops  = 11*n
      nloads  = cprec*0
      nstores = cprec*0

    !MIGHT_BE_BUG: statistics is missing
    case(STAT_PRECTRIGS_BAND_Q1_2D) 
      ! prec_perform_trigs -> PRECS_LOWERLINE_E11_2D in sbblas
      nflops  = 11*n
      nloads  = cprec*0
      nstores = cprec*0

    case(STAT_PRECILU_BAND_Q1_2D) 
      ! one ILU preconditioning step (corresponding routine: if117.f)
      nflops  = (n-1)*9 + 1 + (n-1)*10
      nloads  = cprec*( (n-1)*19 + 3 + (n-1)*20 )
      nstores = cprec*( 2*(n-1) )

    !MIGHT_BE_BUG: statistics is missing
    case(STAT_UMFPACK_2D) 
      ! LU solver from umfpack
      nflops  = 0*n
      nloads  = cprec*0
      nstores = cprec*0

    case(STAT_PROLONGATE_Q1_2D)
      !inner element:     
      !           ______   ______
      !          /      \ /      \
      !        
      !         O________X________O
      !         |\       |       /|  \
      !         | \      |      / |  |
      !         |  \     |     /  |  |
      !  this   |   \    |    /   |  |
      !  side   |    \   |   /    |  |
      !   is    |     \  |  /     |  |
      ! treated |      \ | /      |  |
      ! by the  |        |        |  /
      !  left   X        X        X
      ! adjacent|        |        |
      ! element |      / | \      |  \
      !         |     /  |  \     |  |
      !         |    /   |   \    |  |
      !         |   /    |    \   |  |
      !         |  /     |     \  |  |
      !         | /      |      \ |  |
      !         |/       |       \|  /
      !         O________X________O
      ! This side is treated by the bottom element

      !the boundary element require additional effort, as 2 operations
      !per element occur additionaly. This is NEGLECTED here.
      nflops  = 8 * n ! simple heuristics: fine = coarse*4
      nloads  = cprec*(8*n)
      nstores = cprec*n*4

    case(STAT_RESTRICT_Q1_2D)
      !MIGHT_BE_BUG
!       ioper = 5 * n + 4 * 5 + 4 + 4 * 5 * sqrt(dble(n))
      nflops  = 5 * n * 4 ! simple heuristics: fine = coarse*4
      nloads  = cprec*(5*n*4)
      nstores = cprec*n

    !MIGHT_BE_BUG stat missing
    case(STAT_PROLONGATE_Q2_2D)
      nflops  = 179 * n
      nloads  = 0
      nstores = 0

    !MIGHT_BE_BUG stat missing
    case(STAT_RESTRICT_Q2_2D)
      nflops  = 0
      nloads  = 0
      nstores = 0

    !MIGHT_BE_BUG stat missing
    case(STAT_PROLONGATE_Q2L_2D)
      nflops  = 60 * n
      nloads  = 0
      nstores = 0

    !MIGHT_BE_BUG stat missing
    case(STAT_RESTRICT_Q2L_2D)
      nflops  = 0
      nloads  = 0
      nstores = 0

    case default
      nflops  = 0
      nloads  = 0
      nstores = 0
!      call error_print(ERR_STAT_INVALID_OPCODE, "stat_getcost", &
!                       ERR_NOT_CRITICAL, iarg1=ctype )

    end select

  end subroutine stat_getCost


! ****************************************************************************************


!<subroutine>
  subroutine stat_logOperation_direct(nflops, nloads, nstores)
  !<description>
  ! Adds the given costs to all currently enabled performance counters.
  ! It is recommended to use stat_logOperation() instead if at all possible.
  !</description>
  !<input>
    integer, intent(in) :: nflops
    integer, intent(in) :: nloads
    integer, intent(in) :: nstores
  !</input>
  !<errors>
  ! none
  !</errors>
!</subroutine>

    ! local vars
    type(t_perfCounter),pointer :: rpcPointer

    ! first, all obligatory statistics
    if (btest(stat_RperfCounters(STAT_STATICPERFCOUNTERS)%Ibits, & 
              STAT_PERFCOUNTER_SCARC_COMPLETE)) then
      rpcPointer => stat_RperfCounters(STAT_STATICPERFCOUNTERS) &
                      %Rpc(STAT_PERFCOUNTER_SCARC_COMPLETE)
      rpcPointer%nflops  = rpcPointer%nflops  + nflops
      rpcPointer%nloads  = rpcPointer%nloads  + nloads
      rpcPointer%nstores = rpcPointer%nstores + nstores
    endif

    ! TODO: ifdef this with DISABLE_MULTIDIM
    if (btest(stat_RperfCounters(STAT_STATICPERFCOUNTERS)%Ibits, & 
              STAT_PERFCOUNTER_MDIM_COMPLETE)) then
      rpcPointer => stat_RperfCounters(STAT_STATICPERFCOUNTERS) &
                      %Rpc(STAT_PERFCOUNTER_MDIM_COMPLETE)
      rpcPointer%nflops  = rpcPointer%nflops  + nflops
      rpcPointer%nloads  = rpcPointer%nloads  + nloads
      rpcPointer%nstores = rpcPointer%nstores + nstores
    endif


  end subroutine stat_logOperation_direct


! ****************************************************************************************


!<subroutine>
  subroutine stat_logOperation(ctype, n, cprec, bconst)
  !<description>
  ! Adds the costs associated with the given parameters to all currently enabled
  ! performance counters.
  !</description>
  !<input>
    ! identifier / opcode
    integer, intent(in) :: ctype
    ! vector length
    integer, intent(in) :: n
    ! precision (STAT_SIZEOFFLOAT or STAT_SIZEOFDOUBLE)
    integer, intent(in) :: cprec
    ! constmatmul or varmatmul
    logical, intent(in), optional :: bconst
  !</input>
  !<errors>
  ! none
  !</errors>
!</subroutine>

    integer :: nflops, nloads, nstores

    ! get costs associated with the given operation
    if (present(bconst)) then
      call stat_getCost(ctype, n, cprec, bconst, nflops, nloads, nstores)
    else
      call stat_getCost(ctype, n, cprec, .FALSE., nflops, nloads, nstores)
    endif

    call stat_logOperation_direct(nflops, nloads, nstores)

  end subroutine stat_logOperation


! ****************************************************************************************


!<function>
  function stat_ngetFlops_byID(cperfCounterID) result (nflops)
  !<description>
  ! Returns the number of flops stored in the given static perfcounter.
  ! If the ID represents a dynamic perfcounter, the flops stored in the most 
  ! recently activated one are returned.
  !</description>
  !<input>
    integer, intent(in) :: cperfCounterID
  !</input>
  !<result>
    integer(I64) :: nflops
  !</result>
  !<errors>
  ! ERR_STAT_INVALID_STATICCOUNTER if perfcounter is unknown (critical)
  ! ERR_STAT_INVALID_DYNCOUNTER if perfcounter is unknown (critical)
  !</errors>
!</function>

    ! local vars
    integer :: iarray

    ! figure out which array to work on
    if (cperfCounterID .eq. STAT_PERFCOUNTER_SCARC_DYN) then
      iarray = STAT_SCARCPERFCOUNTERS
    else if (cperfCounterID .eq. STAT_PERFCOUNTER_MULTIDIM_DYN) then
      iarray = STAT_MULTIDIMPERFCOUNTERS
    else 
      iarray = STAT_STATICPERFCOUNTERS  
    endif

    ! return nflops
    if (iarray .eq. STAT_STATICPERFCOUNTERS) then
      nflops = stat_RperfCounters(iarray)%Rpc(cperfCounterID)%nflops
    else
      nflops = stat_RperfCounters(iarray)%Rpc(stat_RperfCounters(iarray)%inext-1)%nflops
    end if

  end function stat_ngetFlops_byID


! ****************************************************************************************


!<function>
  function stat_ngetFlops_byRef(rpc) result (nflops)
  !<description>
  ! Returns the number of flops stored in the given static perfcounter.
  ! If the ID represents a dynamic perfcounter, the flops stored in the most 
  ! recently activated one are returned.
  !</description>
  !<input>
    type(t_perfCounter), intent(in) :: rpc
  !</input>
  !<result>
    integer(I64) :: nflops
  !</result>
  !<errors>
  !  none
  !</errors>
!</function>

    ! return nflops
    nflops = rpc%nflops

  end function stat_ngetFlops_byRef


! ****************************************************************************************


!<function>
  function stat_dgetElapsedCPU(cperfCounterID) result (delapsed)
  !<description>
  ! Returns elapsed time stored in the given perfcounter.
  !</description>
  !<input>
    integer :: cperfCounterID
  !</input>
  !<result>
    real(DP) :: delapsed
  !</result>
  !<errors>
  ! ERR_STAT_INVALID_STATICCOUNTER if perfcounter is unknown (critical)
  ! ERR_STAT_INVALID_DYNCOUNTER if perfcounter is unknown (critical)
  !</errors>
!</function>

    ! local vars
    integer :: iarray

    ! figure out which array to work on
    if (cperfCounterID .eq. STAT_PERFCOUNTER_SCARC_DYN) then
      iarray = STAT_SCARCPERFCOUNTERS
    else if (cperfCounterID .eq. STAT_PERFCOUNTER_MULTIDIM_DYN) then
      iarray = STAT_MULTIDIMPERFCOUNTERS
    else 
      iarray = STAT_STATICPERFCOUNTERS  
    endif

    ! return elapsedtime
    if (iarray .eq. STAT_STATICPERFCOUNTERS) then
      delapsed = stat_RperfCounters(iarray)%Rpc(cperfCounterID)%rtimer%delapsedCPU
    else
      delapsed = stat_RperfCounters(iarray)%Rpc(stat_RperfCounters(iarray)%inext-1) &
                                         %rtimer%delapsedCPU
    end if

  end function stat_dgetElapsedCPU


! ****************************************************************************************


!<function>
  function stat_dgetElapsedReal(cperfCounterID) result (delapsed)
  !<description>
  ! Returns elapsed time stored in the given perfcounter.
  !</description>
  !<input>
    integer :: cperfCounterID
  !</input>
  !<result>
    real(DP) :: delapsed
  !</result>
  !<errors>
  ! ERR_STAT_INVALID_STATICCOUNTER if perfcounter is unknown (critical)
  ! ERR_STAT_INVALID_DYNCOUNTER if perfcounter is unknown (critical)
  !</errors>
!</function>

    ! local vars
    integer :: iarray

    ! figure out which array to work on
    if (cperfCounterID .eq. STAT_PERFCOUNTER_SCARC_DYN) then
      iarray = STAT_SCARCPERFCOUNTERS
    else if (cperfCounterID .eq. STAT_PERFCOUNTER_MULTIDIM_DYN) then
      iarray = STAT_MULTIDIMPERFCOUNTERS
    else 
      iarray = STAT_STATICPERFCOUNTERS  
    endif


    ! return elapsedtime
    if (iarray .eq. STAT_STATICPERFCOUNTERS) then
      delapsed = stat_RperfCounters(iarray)%Rpc(cperfCounterID)%rtimer%delapsedReal
    else
      delapsed = stat_RperfCounters(iarray)%Rpc(stat_RperfCounters(iarray)%inext-1) &
                                         %rtimer%delapsedReal
    end if

  end function stat_dgetElapsedReal


! ****************************************************************************************


!<function>
  function stat_sgetTime_byPCID(cperfCounterID) result (stime)
  ! <description>
  ! Returns a string (length 14) containing both CPU and real elapsed time 
  ! for the given perfcounter.
  ! </description>
  ! <input>
    integer :: cperfCounterID
  ! </input>
  ! <result>
    character(len=14) :: stime
  ! </result>
!</function>

    character(len=6) saux1, saux2

    if (stat_RperfCounters(STAT_STATICPERFCOUNTERS)%Rpc(cperfCounterID)% &
        rtimer%dstartCPU .eq. STAT_TIMER_NOT_RUNNING) then
      saux1 = sys_sdL(stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                      Rpc(cperfCounterID)%rtimer%delapsedCPU, 2)
      saux2 = sys_sdL(stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                      Rpc(cperfCounterID)%rtimer%delapsedReal, 2)
      stime = adjustr(saux1) // " /" // adjustr(saux2)
    else
       stime = "XXXXXXXXXX"
    end if

end function stat_sgetTime_byPCID


! ****************************************************************************************


!<function>
  function stat_sgetTime_byPCRef(rpc) result (stime)
  ! <description>
  ! Returns a string (length 14) containing both CPU and real elapsed time 
  ! for the given perfcounter.
  ! </description>
  ! <input>
    type(t_perfCounter) :: rpc
  ! </input>
  ! <result>
    character(len=14) :: stime
  ! </result>
!</function>

    character(len=6) saux1, saux2

    saux1 = sys_sdL(rpc%rtimer%delapsedCPU, 2)
    saux2 = sys_sdL(rpc%rtimer%delapsedReal, 2)
    stime = adjustr(saux1) // " /" // adjustr(saux2)

end function stat_sgetTime_byPCRef


! ****************************************************************************************


!<function>
  function stat_sgetTimeRatio_byPCRef(rpcP, rpcG) result (stimeRatio)
  ! <description>
!missing
  ! </description>
  ! <input>
    ! missing
    type(t_perfCounter) :: rpcP
    ! missing
    type(t_perfCounter) :: rpcG

  ! </input>
  ! <result>
    character(len=14) :: stimeRatio
  ! </result>
!</function>

    if ((rpcG%rtimer%delapsedCPU .le. 0.0) .or. (rpcG%rtimer%delapsedReal .le. 0.0)) then
      stimeRatio="--/--"
    else
      stimeRatio = adjustr(trim(sys_sdL( &
                     rpcP%rtimer%delapsedCPU / rpcG%rtimer%delapsedCPU * 100.0, 2))) // " % / " // &
                   adjustr(trim(sys_sdL( &
                     rpcP%rtimer%delapsedReal / rpcG%rtimer%delapsedReal * 100.0, 2))) // " %"
    endif


end function stat_sgetTimeRatio_byPCRef


! ****************************************************************************************


!<function>
  function stat_sgetMFlopSec_byID(cperfCounterID) result (smflops)
  ! <description>
  ! Returns a string containing the elapsed floating point operations for the 
  ! given static perfcounter.
  ! </description>
  ! <input>
    integer :: cperfCounterID
  ! </input>
  ! <result>
    character(len=12) :: smflops
  ! </result>
!</function>

    character(len=5) saux1, saux2

    if (stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
        Rpc(cperfCounterID)%rtimer%dstartCPU .eq. STAT_TIMER_NOT_RUNNING) then
      if (stat_RperfCounters(STAT_STATICPERFCOUNTERS)%Rpc(cperfCounterID)% &
          rtimer%delapsedCPU .gt. 0.01) then
        saux1 = sys_si(nint(stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                       Rpc(cperfCounterID)%nflops / &
                       stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                       Rpc(cperfCounterID)%rtimer%delapsedCPU / 1.0E6_DP),5)
      else
        saux1 = "   --"
      endif
      if (stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
          Rpc(cperfCounterID)%rtimer%delapsedReal .gt. 0.01) then
        saux2 = sys_si(nint(stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                       Rpc(cperfCounterID)%nflops / &
                       stat_RperfCounters(STAT_STATICPERFCOUNTERS)% &
                       Rpc(cperfCounterID)%rtimer%delapsedReal / 1.0E6_DP),5)
      else
        saux2 = "   --"
      endif
    else
      saux1 = "   XX"
      saux2 = "   XX"
   endif

   smflops = adjustr(saux1) // " /" // adjustr(saux2)

  end function stat_sgetMFlopSec_byID


! ****************************************************************************************


!<function>
  function stat_sgetMFlopSec_byRef(rpc) result (smflops)
  ! <description>
  ! Returns a string containing the elapsed floating point operations for 
  ! the given perfcounter.
  ! </description>
  ! <input>
    type(t_perfcounter) :: rpc
  ! </input>
  ! <result>
    character(len=12) :: smflops
  ! </result>
!</function>

    character(len=5) saux1, saux2



    if (rpc%rtimer%dstartCPU .eq. STAT_TIMER_NOT_RUNNING) then
      if (rpc%rtimer%delapsedCPU .gt. 0.01) then
        saux1 = sys_si(nint(rpc%nflops / rpc%rtimer%delapsedCPU / 1.0E6_DP),5)
      else
        saux1 = "   --"
      endif
      if (rpc%rtimer%delapsedReal .gt. 0.01) then
        saux2 = sys_si(nint(rpc%nflops / rpc%rtimer%delapsedReal / 1.0E6_DP),5)
      else
        saux2 = "   --"
      endif
    else
      saux1 = "   XX"
      saux2 = "   XX"
   endif

   smflops = adjustr(saux1) // " /" // adjustr(saux2)

  end function stat_sgetMFlopSec_byRef

end module
