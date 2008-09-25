!##############################################################################
!# ****************************************************************************
!# <name> fsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains system routines like time measurement,          
!# string/value conversions and auxiliary routines.
!#
!# On start of the main program, the routine system_init() must be called
!# once to initialise internal values!
!#
!# The following routines can be found here
!#
!#  1.) sys_initClock
!#      -> Initialise the time measurement
!#
!#  2.) sys_doneClock
!#      -> Finalise the time measurement
!#
!#  3.) sys_setClock
!#      -> Set clock for time measurement
!#
!#  4.) sys_setEncompClock
!#      -> Set encompassing clock for time measurement
!#
!#  5.) sys_startClock
!#      -> Start clock for time measurement
!#
!#  6.) sys_stopClock
!#      -> Stop clock for time measurement
!#
!#  7.) sys_stopClockAll
!#      -> Stop all clocks for time measurement
!#
!#  8.) sys_infoClock
!#      -> Output information about time measurement
!#
!#   9.) sys_permute
!#      -> Compute a random permutation of a given sequence
!#
!# 10.) sys_halt
!#      -> Halts the application. Replacement for CALL sys_halt() in F90.
!#         Can be configured how to halt.
!#         E.g. if the global variable sys_haltmode is set to SYS_HALT_THROFPE,
!#         this routine will stop the program by a floating point exception,
!#         which prints the stack trace to the terminal on some compilers.
!#
!# 11.) sys_throwFPE
!#      -> Throw a floating point exception
!#
!# 12.) system_init = system_init_simple /
!#                    system_init_ext
!#      -> Initialise system-wide settings
!#
!# 13.) sys_version
!#      -> Get kernal version number
!#
!# 14.) sys_toupper = sys_toupper_replace /
!#                    sys_toupper_copy
!#      -> Convert a string to uppercase
!#
!# 15.) sys_tolower = sys_tolower_replace /
!#                    sys_tolower_copy
!#      -> Convert a string to lowercase
!#
!# 16.) sys_upcase
!#      -> Convert a string to uppercase, function version
!#
!# 17.) sys_lowcase
!#      -> Convert a string to lowercase, function version
!#
!# 18.) sys_charreplace
!#      -> Replaces characters in a string
!#
!# 19.) sys_getFreeUnit
!#      -> Determine a free file handle for use in an OPEN() command
!#
!# 20.) sys_fileExists
!#      -> Check if file with a given name exists
!#
!# 21.) sys_flush
!#      -> Flush file (if available)
!#
!# 22.) sys_str2Double
!#      -> Convert string to double value
!#
!# 23.) sys_str2Single
!#      -> Convert string to single value
!#
!# </purpose>
!##############################################################################

module fsystem

  implicit none

!<constants>

!<constantblock description="constants for logical values">

  ! logical value 'true'
  integer, parameter :: YES = 0

  ! logical value 'false'
  integer, parameter :: NO = 1

!</constantblock>
  
  integer, parameter :: DP = selected_real_kind(15,307)
  integer, parameter :: SP = selected_real_kind(6,37)

!<constantblock description="Kind values for integers">

  ! kind value for 32Bit integer
  integer, parameter :: I32 = selected_int_kind(8)

  ! kind value for 64Bit integer
  integer, parameter :: I64 = selected_int_kind(10)

!</constantblock>

!<constantblock description="system flags">

  ! constant for a system beep
  character(len=1), parameter :: BEEP = achar(7)

  ! constant for breaking line in a string
  character(len=1), parameter :: NEWLINE = achar(10)

  ! standard length for strings in FEAT
  integer, parameter :: SYS_STRLEN = 256

  ! standard length for name strings in FEAT
  integer, parameter :: SYS_NAMELEN = 32

  ! minimal difference to unity for real values
  real(DP)           :: SYS_EPSREAL = epsilon(1.0_DP)

  ! minimal positive values for real variables
  real(DP)           :: SYS_MINREAL = tiny(1.0_DP)

  ! maximal values for real variables
  real(DP)           :: SYS_MAXREAL = huge(1.0_DP)

  ! maximal values for integer variables
  integer            :: SYS_MAXINT = huge(1)

  ! mathematical constant Pi
  real(DP)           :: SYS_PI
  
  ! internal constant for infinity
  real(DP), parameter :: SYS_INFINITY = huge(1.0_DP)

  ! increment value = 1
  integer, parameter :: INCX = 1

  ! flag for appending data to a file
  integer, parameter :: SYS_APPEND = 0

  ! flag for replacing a file
  integer, parameter :: SYS_REPLACE = 1

!</constantblock>
  
!<constantblock description="system signals">

  integer, parameter :: SIGILL  =  4
  integer, parameter :: SIGTRAP =  5
  integer, parameter :: SIGABRT =  6
  integer, parameter :: SIGEMT  =  7
  integer, parameter :: SIGFPE  =  8
  integer, parameter :: SIGBUS  = 10
  integer, parameter :: SIGSEGV = 11
  
!</constantblock>

!<constantblock description="Constants for the sys_haltmode variable">

  ! Halts the program by the CALL sys_halt() command
  integer, parameter :: SYS_HALT_STOP     = 0

  ! Halts the program by sys_throwFPE. On some compilers, this helps with
  ! debugging as the compiler will print a stack trace to the terminal
  ! that allows tracing back where an error came from.
  integer, parameter :: SYS_HALT_THROWFPE = 1
  
!</constantblock>

!</constants>

!************************************************************************

!<types>

!<typeblock>

  ! Emulation of an array of pointers to double precision vectors
  type t_realPointer
    real(DP), dimension(:), pointer :: ptr
  end type

!</typeblock>

!<typeblock>

  ! Emulation of an array of pointers to 2D double precision vectors
  type t_realw2DPointer
    real(DP), dimension(:,:), pointer :: ptr
  end type

!</typeblock>

!<typeblock>

  ! Emulation of an array of pointers to single precision vectors
  type t_singlePointer
    real(SP), dimension(:), pointer :: ptr
  end type

!</typeblock>
  
!<typeblock>

  ! Emulation of an array of pointers to 2D single precision vectors
  type t_single2DPointer
    real(SP), dimension(:), pointer :: ptr
  end type

!</typeblock>

!<typeblock>

  ! Emulation of an array of pointers to integer vectors
  type t_intPointer
    integer(I32), dimension(:), pointer :: ptr
  end type
  
!</typeblock>
  
!<typeblock>

  ! Emulation of an array of pointers to 2D integer vectors
  type t_int2DPointer
    integer(I32), dimension(:,:), pointer :: ptr
  end type
  
!</typeblock>
  
!<typeblock>

  ! Global project settings and system information
  type t_sysconfig

    ! project id
    character(LEN=SYS_STRLEN) :: sprojectID    = ''

    ! project directory
    character(LEN=SYS_STRLEN) :: sprojectDir   = ''

    ! starting time of this project
    integer                   :: iprojectStart = 0

    ! starting time of this project (long time runs)
    integer, dimension(8)     :: iprojectStartLong = 0
    
  end type

!</typeblock>

!<typeblock>

  ! Global time measurement
  type t_clock

    ! The name of the clock
    character(LEN=SYS_NAMELEN) :: sname = ''

    ! Starting time
    integer  :: istart       = 0

    ! Stopping time
    integer  :: istop        = 0

    ! Total time elapsed
    integer :: icount        = 0

    ! Total overflow of system clock
    integer :: ioverflow     = 0

    ! Activation flag
    logical  :: bactive      = .false.

    ! Number of encompassing clock
    integer  :: iencompClock = 0

  end type t_clock

!</typeblock>

!</types>

!************************************************************************

!<globals>
  ! Global time measurement structure
  type(t_clock), dimension(:), allocatable, target, save :: rclock

  ! maximal measurable time span in seconds (system-dependend)
  real(DP) :: sys_dtimeMax

  ! global system configuration
  type (t_sysconfig), target, save :: sys_sysconfig
  
  ! Halt mode. This variable defines the way, sys_halt halts the program.
  ! One of the SYS_HALT_xxxx constants.
  integer :: sys_haltmode = SYS_HALT_STOP
!</globals>

!************************************************************************

  interface system_init
    module procedure system_init_simple
    module procedure system_init_ext
  end interface

  interface sys_toupper
    module procedure sys_toupper_replace
    module procedure sys_toupper_copy
  end interface

  interface sys_tolower
    module procedure sys_tolower_replace
    module procedure sys_tolower_copy
  end interface

contains

!************************************************************************

!<subroutine>

  subroutine sys_initClock(nclock)

!<description>

    ! This routine initializes the time measurement.

!</description>

!<input>

    ! Number of clocks.
    integer, intent(IN) :: nclock

!</input>

!</subroutine>

    if (allocated(rclock)) then
      print *, "sys_initClock: Time measurement is already initialised!"
      call sys_halt()
    end if
    
    allocate(rclock(nclock))
    
  end subroutine sys_initClock

!************************************************************************

!<subroutine>

  subroutine sys_doneClock

!<description>

    ! This routine releases the time measurement.

!</description>

!</subroutine>

    if (.not.allocated(rclock)) then
      print *, "sys_doneClock: Time measurement is not initialised!"
      call sys_halt()
    end if
    
    deallocate(rclock)
    
  end subroutine sys_doneClock


!************************************************************************

!<subroutine>

  subroutine sys_setClock(iclock, sname)

!<description>

  ! This routine sets a clock for time measurement.

!</description>

!<input>

    ! Number of the clock
    integer, intent(IN)          :: iclock

    ! Name of the clock
    character(LEN=*), intent(IN) :: sname

!</input>

!</subroutine>

    if (.not.allocated(rclock)) then
      print *, "sys_setClock: Time measurement is not initialised!"
      call sys_halt()
    end if
    
    if (iclock .gt. size(rclock)) then
      print *, "sys_setClock: Clock number exceeds maximum number of clocks!"
      call sys_halt()
    end if
    
    rclock(iclock)%sname        = sname
    rclock(iclock)%istart       = 0
    rclock(iclock)%istop        = 0
    rclock(iclock)%icount       = 0
    rclock(iclock)%ioverflow    = 0
    rclock(iclock)%bactive      = .false.
    rclock(iclock)%iencompClock = 0

  end subroutine sys_setClock

!************************************************************************

!<subroutine>

  subroutine sys_setEncompClock(iclock, iencompClock)

!<description>

  ! This routine sets an encompassing clock for a given clock.

!</description>

!<input>

    ! Number of the clock
    integer, intent(IN) :: iclock

    ! Number of the encompassing clock
    integer, intent(IN) :: iencompClock

!</input>
    
!</subroutine>

    if (.not.allocated(rclock)) then
      print *, "sys_setEncompClock: Time measurement is not initialised!"
      call sys_halt()
    end if
    
    if (max(iclock,iencompClock) .gt. size(rclock)) then
      print *, "sys_setEncompClock: Clock number(s) exceed maximum number of clocks!"
      call sys_halt()
    end if
    
    rclock(iclock)%iencompClock = iencompClock
    
  end subroutine sys_setEncompClock

!************************************************************************

!<subroutine>

  subroutine sys_startClock(iclock)

!<description>

    ! This routine starts the given clock.

!</description>

!<input>

    ! Number of the clock
    integer, intent(IN) :: iclock

!</input>

!</subroutine>

    if (.not.allocated(rclock)) then
      print *, "sys_startClock: Time measurement is not initialised!"
      call sys_halt()
    end if
    
    if (iclock .gt. size(rclock)) then
      print *, "sys_startClock: Clock number exceeds maximum number of clocks!"
      call sys_halt()
    end if
    
    ! Start time measurement
    rclock(iclock)%bactive = .true.
    call system_clock(rclock(iclock)%istart)
  end subroutine sys_startClock

!************************************************************************

!<subroutine>

  subroutine sys_stopClock(iclock)

!<description>

    ! This routine stops the given clock.

!</description>

!<input>

    ! Number of the clock
    integer, intent(IN) :: iclock

!</input>

!</subroutine>

    ! local variables
    integer :: icount,irate,icountmax

    if (.not.allocated(rclock)) then
      print *, "sys_stopClock: Time measurement is not initialised!"
      call sys_halt()
    end if
    
    if (iclock .gt. size(rclock)) then
      print *, "sys_stopClock: Clock number exceeds maximum number of clocks!"
      call sys_halt()
    end if
    
    if (.not.rclock(iclock)%bactive) then
      print *, "sys_stopClock: Clock has not been started, skipping!"
      return
    end if

    ! Stop time measurement
    rclock(iclock)%bactive = .false.
    call system_clock(rclock(iclock)%istop, irate, icountmax)
    
    ! Calculate elapsed time. Note that the stopping time may be smaller
    ! than the starting time. In this case the maximum number of counts
    ! has been reached and must be considered.
    if (rclock(iclock)%istart  .gt. rclock(iclock)%istop) then
      icount = icountmax-rclock(iclock)%istart
      icount = icount+rclock(iclock)%istop
      if (rclock(iclock)%icount .gt. icountmax-icount) then
        rclock(iclock)%icount    = rclock(iclock)%icount-icountmax
        rclock(iclock)%icount    = rclock(iclock)%icount+icount
        rclock(iclock)%ioverflow = rclock(iclock)%ioverflow+1
      else
        rclock(iclock)%icount =  rclock(iclock)%icount+icount
      end if
    else
      icount = rclock(iclock)%istop-rclock(iclock)%istart
      if (rclock(iclock)%icount .gt. icountmax-icount) then
        rclock(iclock)%icount    = rclock(iclock)%icount-icountmax
        rclock(iclock)%icount    = rclock(iclock)%icount+icount
        rclock(iclock)%ioverflow = rclock(iclock)%ioverflow+1
      else
        rclock(iclock)%icount = rclock(iclock)%icount+icount
      end if
    end if
  end subroutine sys_stopClock

!************************************************************************

!<subroutine>

  subroutine sys_stopClockAll

!<description>

    ! This routine stops all clocks.

!</description>

!</subroutine>

    ! local variables
    integer :: iclock

    if (.not.allocated(rclock)) then
      print *, "sys_stopClockAll: Time measurement is not initialised!"
      call sys_halt()
    end if

    do iclock = 1, size(rclock)
      call sys_stopClock(iclock)
    end do

  end subroutine sys_stopClockAll

!************************************************************************

!<subroutine>

  subroutine sys_infoClock

!<description>

    ! This routine prints information about the clock

!</description>

!</subroutine>

    ! local variables
    type(t_clock), dimension(:), pointer :: rclockTmp
    real(DP) :: dtotaltime,dtime,dtimeself
    integer  :: iclock,icount,irate,icountmax
    integer  :: ndays,nsecs,nsec1,nsec2
    integer, dimension(8) :: istopLong

    write(*,FMT='(A)') 'Time measurement:'
    write(*,FMT='(A)') '-----------------'

    ! Compute total cpu time
    call system_clock(icount,irate,icountmax)
    call date_and_time(values=istopLong)

    ! Compute number of days elapsed
    ndays = calender_to_julian(istopLong(1), istopLong(2), istopLong(3))-&
        calender_to_julian(sys_sysconfig%iprojectStartLong(1),&
                           sys_sysconfig%iprojectStartLong(2),&
                           sys_sysconfig%iprojectStartLong(3))

    ! Compute number of seconds elapsed
    nsec1 = istopLong(5)*3600+istopLong(6)*60+istopLong(7)
    nsec2 = sys_sysconfig%iprojectStartLong(5)*3600+&
            sys_sysconfig%iprojectStartLong(6)*60+&
            sys_sysconfig%iprojectStartLong(7)

    if (nsec2 .lt. nsec1) then
      nsecs = 86400-(nsec2-nsec1)
    else
      nsecs = nsec2-nsec1
    end if

    ! Compute total time from system clock
    if (icount .lt. sys_sysconfig%iprojectStart) then
      dtime = real(icountmax+icount-sys_sysconfig%iprojectStart, DP)/real(irate, DP)
    else
      dtime = real(icount-sys_sysconfig%iprojectStart, DP)/real(irate, DP)
    end if

    ! Ok, now we have two times. Check which measurement can be used
    dtotaltime = real(icountmax, DP)/real(irate, DP)

    if (ndays .le. int(dtotaltime/86400._DP)) then
      dtotaltime = dtime
      write(*,FMT='(A)') 'Simulation ran '//trim(sys_sdEL(dtime, 8))//' seconds'
    else
      dtotaltime = real(ndays*86400+nsecs, DP)
      write(*,FMT='(A)') 'Simulation ran '//trim(sys_siL(ndays, 15))//' days and '//&
                                            trim(sys_siL(nsecs, 15))//' seconds'
    end if

    if (.not.allocated(rclock)) then
      print *, "sys_infoClock: Time measurement is not initialised!"
      call sys_halt()
    end if

    ! Make a copy of the system clock
    allocate(rclockTmp(size(rclock)))
    rclockTmp = rclock

    ! Subtract the individual times from the encompassing clocks
    do iclock = 1, size(rclockTmp)
      if (rclockTmp(iclock)%iencompClock .ne. 0) then
        rclockTmp(rclockTmp(iclock)%iencompClock)%icount =&
            rclockTmp(rclockTmp(iclock)%iencompClock)%icount-rclockTmp(iclock)%icount
        if (rclockTmp(rclockTmp(iclock)%iencompClock)%icount .lt. 0) then
          rclockTmp(rclockTmp(iclock)%iencompClock)%icount =&
              rclockTmp(rclockTmp(iclock)%iencompClock)%icount+icountmax
          rclockTmp(rclockTmp(iclock)%iencompClock)%ioverflow =&
              rclockTmp(rclockTmp(iclock)%iencompClock)%ioverflow-1  
        end if
      end if
    end do
    
    ! Check consistency
    do iclock = 1, size(rclockTmp)
      if (rclockTmp(iclock)%icount .lt. 0 .or.&
          rclockTmp(iclock)%ioverflow .lt. 0) then
        print *, "sys_infoClock: Time measurement is incorrect!"
      end if
    end do

    write(*,*)
    write(*,FMT='(A,T35,A,T60,A)') 'Clock name', 'elapsed time (%)', 'self time (%)'
    write(*,FMT='(85("="))')

    ! Print out all clocks
    do iclock = 1, size(rclock)
      
      ! Compute time in seconds
      dtime = rclock(iclock)%ioverflow*real(icountmax, DP)/real(irate, DP)+&
              real(rclock(iclock)%icount, DP)/real(irate, DP)
      dtimeself = rclockTmp(iclock)%ioverflow*real(icountmax, DP)/real(irate, DP)+&
                  real(rclockTmp(iclock)%icount, DP)/real(irate, DP)

      write(*,FMT='(A,T35,A,T60,A)')&
          trim(adjustl(rclock(iclock)%sname)),&
          trim(sys_sdEL(dtime, 8))//' ('//&
          trim(sys_sdL(100._DP/dtotaltime*dtime, 2))//')',&
          trim(sys_sdEL(dtimeself, 8))//' ('//&
          trim(sys_sdL(100._DP/dtotaltime*dtimeself,2))//')'
    end do

    ! Free unused memory
    deallocate(rclockTmp)

  contains

    function calender_to_julian(year, month, day) result(ivalue)
      integer, intent(IN) :: year
      integer, intent(IN) :: month
      integer, intent(IN) :: day

      integer             :: ivalue

      ivalue = day-32075+&
               1461*(year+4800+(month-14)/12)/4+&
               367*(month-2-((month-14)/12)*12)/12-&
               3*((year+4900+(month-14)/12)/100)/4
    end function calender_to_julian
  end subroutine sys_infoClock

!************************************************************************

!<subroutine>

  subroutine sys_permute(k,Idata)

!<description>
    ! This routine computes the permutation of the initial set Idata
    ! that corresponds to the factoriodic number k.
!</description>

!<input>
    ! factoriodic number
    integer(I32), intent(IN) :: k
!</input>

!<inputoutput>
    ! initial and permuted set
    integer(I32), dimension(:), intent(INOUT) :: Idata
!</inputoutput>
!</subroutine>
  
    ! local variables
    integer(I32) :: factorial,i,j,iswap
    
    do j=2,size(Idata)

      if (j-1 .lt. k) then
        factorial = 0._I64
      elseif(j-1 .eq. k) then
        factorial = 1._I64
      else
        factorial=k+1
        do i=k+2,j-1
          factorial=factorial*i
        end do
      end if
      
      i = j-mod(factorial,j)
      iswap    = Idata(j)
      Idata(j) = Idata(i)
      IData(i) = iswap
    end do
  end subroutine sys_permute

!************************************************************************

!<subroutine>

  subroutine sys_halt()
  
!<description>
    ! This routine halts the application like the CALL sys_halt() command in
    ! Fortran 90. The routine can be configured how to halt the application.
    ! For this purpose, the main program can set the global variable
    ! sys_haltmode to one of the SYS_HALT_xxxx constants.
!</description>
    
!</subroutine>

    select case (sys_haltmode)
    case (SYS_HALT_STOP)
      stop
    case (SYS_HALT_THROWFPE)
      call sys_throwFPE()
    end select

  end subroutine sys_halt

!************************************************************************

!<subroutine>

  pure subroutine sys_throwFPE()
  
!<description>
    ! This routine throws a floating point exception for debugging
    ! purposes to prevent the debugger to exit the program.
!</description>
    
!</subroutine>
    
    ! local variables
    integer :: i1,i2

    i1=1
    i2=0
    i1=i1/i2
    
  end subroutine sys_throwFPE

!************************************************************************

!<subroutine>

  subroutine system_init_simple()

!<description>
    ! This subroutine initialises internal data structures 
    ! with standard values.
!</description>

!</subroutine>

    call system_init_ext("","")
    
  end subroutine system_init_simple

!************************************************************************

!<subroutine>

  subroutine system_init_ext(sprojectID,sprojectDir)

!<description>
    ! Extended initialisation.
    ! This subroutine initialises internal data structures.
    ! The value of the project name and directory are set according to the
    ! parameters.
!</description>

!<input>
    ! An ID string of the project.
    character(LEN=*), intent(IN) :: sprojectID
    
    ! The directory of the project. "" means 'current directory'.
    character(LEN=*), intent(IN) :: sprojectDir
!</input>

!</subroutine>

    ! local variables
    integer :: icount ! current system time
    integer :: irate  ! approx. number of system clock ticks per second
    integer :: icmax  ! largest possible value of icount

    ! system_clock is not a FEAT, but a basic FORTRAN 90 routine
    call system_clock(icount,irate,icmax)

    ! use data_and_time to measure long time runs
    call date_and_time(values=sys_sysconfig%iprojectStartLong)

    ! maximal measurable time span in seconds (system-dependend)
    sys_dtimeMax = real(icmax,DP)/real(irate,DP)

    ! Initialise the global sysconfig structure
    sys_sysconfig%sprojectID    = sprojectID
    sys_sysconfig%sprojectDir   = sprojectDir
    sys_sysconfig%iprojectStart = icount

    ! Set value of Pi = 3.14..
    SYS_PI=asin(1.0_DP)*2.0_DP
    
  end subroutine system_init_ext

!************************************************************************************


!<subroutine>

  subroutine sys_version(ifeastVersionHigh, ifeastVersionMiddle, ifeastVersionLow, &
                         sreldate)

!<description>
    ! This subroutine returns the library version information.
!</description>

!<output>

    ! high version number
    integer :: ifeastVersionHigh
    
    ! middle version number
    integer :: ifeastVersionMiddle
    
    ! low version number
    integer :: ifeastVersionLow
    
    ! release date
    character(LEN=*) :: sreldate

!</output>

!</subroutine>

    ifeastVersionHigh=0
    ifeastVersionMiddle=0
    ifeastVersionLow=1
    
    sreldate="01.01.2007 RC0"
    
  end subroutine sys_version

!************************************************************************

!<subroutine>
  
  subroutine sys_toupper_replace (str)

!<description>
    ! Convert a string to upper case. 
    ! The given string is replaced by its uppercase version.
!</description>

!<inputoutput>
  
    ! The string that is to make uppercase
    character(LEN=*), intent(INOUT) :: str
  
!</inputoutput>
  
!</subroutine>
  
    ! local variables
    integer, parameter :: up2low = iachar("a") - iachar("A")
    integer :: i
    character    :: c
    
    do i=1,len(str)
      c = str(i:i)
      if ((c .ge. "a") .and. (c .le. "z")) then
        str(i:i) = achar (iachar(c) - up2low)
      end if
    end do
    
  end subroutine sys_toupper_replace

!************************************************************************

!<subroutine>
  
  subroutine sys_toupper_copy (str,strUpper) 

!<description>
    ! Convert a string to upper case.
!</description>

!<input>
  
    ! The string that is to make uppercase
    character(LEN=*), intent(IN) :: str

!</input>

!<output>

    ! Uppercase version of the given string
    character(LEN=*), intent(OUT) :: strUpper
  
!</output>
  
!</subroutine>
  
    ! local variables
    integer, parameter :: up2low = iachar("a") - iachar("A")
    integer :: i
    character    :: c
    
    if (len(str) > len(strUpper)) then
      print *, "sys_toupper_copy: target string is too short"
      call sys_halt()
    end if
    
    ! Initialize string
    strUpper = ''
    
    do i=1,len(str)
      c = str(i:i)
      if ((c .ge. "a") .and. (c .le. "z")) then
        strUpper(i:i) = achar (iachar(c) - up2low)
      else
        strUpper(i:i) = c
      end if
    end do
    
  end subroutine sys_toupper_copy

!************************************************************************

!<subroutine>
  
  subroutine sys_tolower_replace (str) 

!<description>
    ! Convert a string to lower case.
    ! The given string is replaced by its lowercase version.
!</description>

!<inputoutput>
  
    ! The string that is to make lowercase
    character(LEN=*), intent(INOUT) :: str

!</inputoutput>
  
!</subroutine>
  
    ! local variables
    integer, parameter :: up2low = iachar("a") - iachar("A")
    integer :: i
    character    :: c
    
    do i=1,len(str)
      c = str(i:i)
      if ((c .ge. "A") .and. (c .le. "Z")) then
        str(i:i) = achar (iachar(c) + up2low)
      end if
    end do
    
  end subroutine sys_tolower_replace
  
!************************************************************************

!<subroutine>
  
  subroutine sys_tolower_copy (str,strLower) 

!<description>
    ! Convert a string to lower case.
!</description>

!<input>
  
    ! The string that is to make lowercase
    character(LEN=*), intent(IN) :: str

!</input>

!<output>

    ! Lowercase version of the given string
    character(LEN=*), intent(OUT) :: strLower
  
!</output>
  
!</subroutine>
  
    ! local variables
    
    integer, parameter :: up2low = iachar("a") - iachar("A")
    integer :: i
    character    :: c
    
    if (len(str) > len(strLower)) then
      print *, "sys_tolower_copy: target string is too short"
      call sys_halt()
    end if
    
    ! Initialize string
    strLower = ''
    
    do i=1,len(str)
      c = str(i:i)
      if ((c .ge. "A") .and. (c .le. "Z")) then
        strLower(i:i) = achar (iachar(c) + up2low)
      else
        strLower(i:i) = c
      end if
    end do
    
  end subroutine sys_tolower_copy
  
!******************************************************************************

!<function>

  pure function sys_upcase(sinput) result(soutput)
  
!<description>
    ! This routine converts a given string to its uppercase version.
!</description>

!<input>

    ! input string
    character(len=*), intent(in) :: sinput

!</input>

!<output>

    ! output string
    character(len=len(sinput)) :: soutput

!</output>
!</function>

    ! index variable
    integer :: i
    
    soutput = " "   ! initialise string
    do I = 1,len(sinput)
      if(sinput(i:i) .ge. "a" .and. sinput(i:i) .le. "z") then
        soutput(i:i) = achar(iachar(sinput(i:i)) - 32)
      else
        soutput(i:i) = sinput(i:i)
      end if
    end do
    
  end function sys_upcase

!******************************************************************************

!<function>

  pure function sys_lowcase(sinput) result(soutput)
  
!<description>
    ! This routine converts a given string to its uppercase version.
!</description>

!<input>

    ! input string
    character(len=*), intent(in) :: sinput

!</input>

!<output>

    ! output string
    character(len=len(sinput)) :: soutput

!</output>
!</function>

    ! index variable
    integer :: i
    
    soutput = " "   ! initialise string
    do i = 1,len(sinput)
      if(sinput(i:i) .ge. "A" .and. sinput(i:i) .le. "Z") then
        soutput(i:i) = achar(iachar(sinput(i:i)) + 32)
      else
        soutput(i:i) = sinput(i:i)
      end if
    end do
    
  end function sys_lowcase

!******************************************************************************

!<function>

  pure function sys_charreplace(sinput,scharsource,schardest) result(soutput)
  
!<description>
    ! Replaces all characers scharsource in sinput by schardest.
    ! Case sensitive.
!</description>

!<input>
    ! input string
    character(LEN=*), intent(IN) :: sinput
    
    ! Character to be searched for.
    character, intent(IN) :: scharsource
    
    ! Detinatiion character, all scarsource characters in sinput should be
    ! replaced by.
    character, intent(IN) :: schardest
!</input>

!<output>
    ! output string
    character(LEN=len(sinput)) :: soutput
!</output>
!</function>

    !index variable
    integer :: i
    
    soutput = " "   !initialise string
    do I = 1,len(sinput)
      if(sinput(i:i) .eq. scharsource) then
        soutput(i:i) = schardest
      else
        soutput(i:i) = sinput(i:i)
      end if
    end do
    
  end function sys_charreplace

!************************************************************************

!<function>

  integer function sys_getFreeUnit()

!<description>
    !This routine tries to find a free unit (for file input/output). If a free unit is
    !found, it is returned, otherwise -1 is returned.
!</description>

!<result>
    !number of free unit (-1 if no free unit available)
!</result>

!</function>

    logical :: bexists, bopened!flags indicating errors
    integer :: itry !free unit candidate
    
    sys_getFreeUnit = -1
    do itry = 20,10000
      !does unit exist?
      inquire(unit=itry, exist=bexists)
      if (bexists) then
        !is unit already opened?
        inquire(unit=itry, opened=bopened)
        if (.not. bopened) then
          !free unit found
          sys_getFreeUnit = itry
          !exit do-loop
          exit
        endif
      endif
    enddo
    if (sys_getFreeUnit .eq. -1) then
      write (6,*) "*** WARNING! No free unit between 1 and 10000 found! ***"
    endif
    
  end function sys_getFreeUnit


!************************************************************************

!<function>

  logical function sys_fileExists(iunit,sname)

!<description>
    !This function checks if there is a file connected to unit iunit, which we
    !can access for reading.
!</description>

!<input>
    
    !unit the file shall be attached to
    integer :: iunit
    
    !name of the file to look at
    character (len=*):: sname
    
!</input>

!<result>
    !  .TRUE. if the file is accessable for reading,
    !  .FALSE. otherwise
!</result>
!</function>

    integer :: iostat !status variable for opening procedure
    
    open(iunit,FILE=sname,IOSTAT=iostat,STATUS='OLD',ACTION='READ')
    sys_fileExists=(iostat .eq. 0)
    close(iunit)

  end function sys_fileExists


!************************************************************************

!<subroutine>
  subroutine sys_flush(iunit)
    
!<description>
    ! This routine flushes the buffers associated with an open output unit.
    ! This normally happens when the file is closed or the program ends, 
    ! but this routine ensures the buffers are flushed before any other 
    ! processing occurs.
!</description>

!<input>

    !unit connected to the file to write to
    integer :: iunit

!</input>
!</subroutine>

!#ifdef HAS_FLUSH
!    call flush(iunit)
!#endif

  end subroutine sys_flush

!************************************************************************

!<function>

  function sys_str2Double(svalue,sformat) result(dvalue)

!<description>
    ! This routine converts a given string that provides a valid
    ! IEEE 745 representation of a real number into a double value.
!</description>

!<input>

    ! string containing the real number
    character(LEN=*), intent(IN) :: svalue

    ! format description to use for conversion
    character(LEN=*), intent(IN) :: sformat

!</input>

!<result>

    ! double precision value
    real(DP) :: dvalue

!</result>
!</function>
    
    ! local variables
    character(LEN=len(svalue)+3) :: svalueTemp
    integer :: ipos

    ! Check if string contains a 'dot'
    if (scan(svalue,'.') .ne. 0) then

      ! Read original string
      read(svalue,sformat) dvalue
    else
      ! Check if string is given in scientific notation
      ipos = scan(svalue,'dDeE')

      if (ipos .eq. 0) then
        ! Append '.E0' to convert string into scientific notation
        svalueTemp = trim(svalue)//'.E0'
       
      elseif (ipos .eq. 1) then
        ! Prepend '1.' to convert string into scientific notation
        svalueTemp = '1.'//adjustl(svalue)
      else
        ! Insert '.' to convert string into scientific notation
        svalueTemp = svalue(1:ipos-1)//'.'//svalue(ipos:)
      end if

      ! Read modified string
      read(svalueTemp,sformat) dvalue
    end if
    
  end function sys_str2Double

!************************************************************************

!<function>

  function sys_str2Single(svalue,sformat) result(fvalue)

!<description>
    ! This routine converts a given string that provides a valid
    ! IEEE 745 representation of a real number into a single value.
!</description>

!<input>

    ! string containing the real number
    character(LEN=*), intent(IN) :: svalue

    ! format description to use for conversion
    character(LEN=*), intent(IN) :: sformat

!</input>

!<result>

    ! single precision value
    real(SP) :: fvalue

!</result>
!</function>
    
    ! local variables
    character(LEN=len(svalue)+3) :: svalueTemp
    integer :: ipos

    ! Check if string contains a 'dot'
    if (scan(svalue,'.') .ne. 0) then

      ! Read original string
      read(svalue,sformat) fvalue
    else
      ! Check if string is given in scientific notation
      ipos = scan(svalue,'dDeE')

      if (ipos .eq. 0) then
        ! Append '.E0' to convert string into scientific notation
        svalueTemp = trim(svalue)//'.E0'
       
      elseif (ipos .eq. 1) then
        ! Prepend '1.' to convert string into scientific notation
        svalueTemp = '1.'//adjustl(svalue)
      else
        ! Insert '.' to convert string into scientific notation
        svalueTemp = svalue(1:ipos-1)//'.'//svalue(ipos:)
      end if

      ! Read modified string
      read(svalueTemp,sformat) fvalue
    end if
    
  end function sys_str2Single


!************************************************************************
! Main conversion routines:
!
! sys_sl  :  logical   => string
! First parameter: logical value to convert,

! sys_sd  :  real      => string,
! sys_sdE :  real      => string (scientific notation)
! First parameter: real value to convert,
! Second paramter: number of decimals places

! sys_si  :  int       => string
! sys_si0 :  int       => string (filled with zeros)
! sys_sli :  long int  => string
! sys_sli0:  long int  => string (filled with zeros)
! First parameter: integer value to convert,
! Second paramter: number of total digits (filled with white spaces or zeros)

! All routines exist also in a *L version which do basically the same,
! but return a left-adjusted string (fixed length of 32 characters)
!************************************************************************

!<function>
  
  character (len=32) function sys_sl(lvalue) result(soutput)

!<description>
    ! This routine converts a logical value to a string.
!</description>

!<input>

    ! value to be converted
    logical, intent(IN) :: lvalue
!</input>
!</function>

    if (lvalue) then
      soutput = 'true'
    else
      soutput = 'false'
    end if
  end function sys_sl

!<function>

  character (len=32) function sys_sd(dvalue, idigits) result(soutput)

!<description>
    ! This routine converts a double value to a string with idigits
    ! decimal places.
!</description>

!<result>
    ! String representation of the value, filled with white spaces.
    ! At most 32 characters supported.
!</result>

!<input>

    ! value to be converted
    real(DP), intent(in) :: dvalue

    !number of decimals
    integer, intent(in)  :: idigits
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux

    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_sd! ***"
      write(saux, '(i2)') 16
    else
      write(saux, '(i2)') idigits
    endif
    
    sformat = "(f32." // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) dvalue
  end function sys_sd

!************************************************************************

!<function>

  character (len=32) function sys_sdP(dvalue, ipositions, idigits) result(soutput)

!<description>
    ! This routine converts a double value to a string with length
    ! iposition and idigits decimal places.
!</description>

!<result>
    ! String representation of the value, filled with white spaces.
    ! At most 32 characters supported.
!</result>

!<input>

    ! value to be converted
    real(DP), intent(in) :: dvalue

    ! number of positions in the string
    integer, intent(in)  :: ipositions

    ! number of decimals
    integer, intent(in)  :: idigits
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux,saux2

    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_sdP! ***"
      write(saux, '(i2)') 16
    else
      write(saux, '(i2)') idigits
    endif

    if (idigits .gt. 32) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_sdP! ***"
      write(saux2, '(i2)') 32
    else
      write(saux2, '(i2)') ipositions
    endif

    sformat = "(f"//trim(saux2)//"." // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) dvalue
  end function sys_sdP


!************************************************************************


!<function>

  character (len=24) function sys_sdE(dvalue, idigits) result(soutput)

!<description>
    ! This routine converts a double value to a string with idigits
    ! decimal places in scientific notation.
!</description>

!<result>
    ! String representation of the value, filled with white spaces.
    ! At most 24 characters supported.
!</result>

!<input>

    ! value to be converted
    real(DP), intent(in) :: dvalue

    !number of decimals
    integer, intent(in)  :: idigits
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux

    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_sdE! ***"
      write(saux, '(i2)') 16
    else
      write(saux, '(i2)') idigits
    endif

    sformat = "(es24." // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) dvalue
  end function sys_sdE


!************************************************************************

!<function>

  character (len=32) function sys_sdEP(dvalue, ipositions, idigits) result(soutput)

!<description>
    ! This routine converts a double value to a string with length
    ! iposition and idigits decimal places in scientific notation.
!</description>

!<result>
    ! String representation of the value, filled with white spaces.
    ! At most 32 characters supported.
!</result>

!<input>

    ! value to be converted
    real(DP), intent(in) :: dvalue

    ! number of positions in the string
    integer, intent(in)  :: ipositions

    ! number of decimals
    integer, intent(in)  :: idigits
    
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux,saux2

    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_sdEP! ***"
      write(saux, '(i2)') 16
    else
      write(saux, '(i2)') idigits
    endif

    if (idigits .gt. 24) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_sdEP! ***"
      write(saux2, '(i2)') 24
    else
      write(saux2, '(i2)') ipositions
    endif

    sformat = "(es"//trim(saux2)//"." // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) dvalue
  end function sys_sdEP

!************************************************************************


!<function>

  character (len=32) function sys_si(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts an integer value to a string of length idigits.
!</description>

!<result>
    ! String representation of the value, filled with white spaces.
    ! At most 32 characters supported.
!</result>

!<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux
    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_si! ***"
      write(saux, '(i2)') 16
    else if (idigits .lt. 10) then
      write(saux, '(i1)') idigits
    else
      write(saux, '(i2)') idigits
    endif

    sformat = "(i" // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) ivalue

  end function sys_si


!************************************************************************


!<function>

  character (len=32) function sys_si0(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts an integer value to a string of length idigits.
!</description>

!<result>
    ! String representation of the value, filled with zeros.
    ! At most 32 characters supported.
!</result>

!<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux

    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_si0! ***"
      write(saux, '(i2)') 16
    else
      write(saux, '(i2)') idigits
    endif

    sformat = "(i" // trim(saux) // "." // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) ivalue
  end function sys_si0


!************************************************************************


!<function>

  character (len=32) function sys_sli(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts a long integer value to a string of length idigits.
!</description>

!<result>
    ! String representation of the value, filled with white spaces.
    ! At most 32 characters supported.
!</result>

!<input>

    ! value to be converted
    integer(I64), intent(in) :: ivalue

    !number of decimals
    integer, intent(in)      :: idigits
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux

    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_sli! ***"
      write(saux, '(i2)') 16
    else
      write(saux, '(i2)') idigits
    endif

    sformat = "(i" // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) ivalue
  end function sys_sli


!************************************************************************


!<function>

  character (len=32) function sys_sli0(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts a long integer value to a string of length idigits.
!</description>

!<result>
    ! String representation of the value, filled with zeros.
    ! At most 32 characters supported.
!</result>

!<input>

    ! value to be converted
    integer(I64), intent(in) :: ivalue

    !number of decimals
    integer, intent(in)      :: idigits
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux

    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      write(6, *) "*** WARNING! Too many decimal places requested in sys_sli0! ***"
      write(saux, '(i2)') 16
    else
      write(saux, '(i2)') idigits
    endif

    sformat = "(i" // trim(saux) // "." // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) ivalue
  end function sys_sli0


!************************************************************************


!************************************************************************
! Left-adjusted versions of the main conversion routines,
! just add capital L to function name
!************************************************************************


!<function>

  character (len=32) function sys_sdL(dvalue, idigits) result(soutput)

!<description>
    ! This routine converts a double value to a string with idigits
    ! decimal places.
!</description>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of 32 characters
!</result>

!<input>

    ! value to be converted
    real(DP), intent(in) :: dvalue

    !number of decimals
    integer, intent(in)  :: idigits
!</input>
!</function>

    soutput = adjustl(sys_sd(dvalue, idigits))
  end function sys_sdL


!************************************************************************


!<function>

  character (len=32) function sys_sdEL(dvalue, idigits) result(soutput)

!<description>
    ! This routine converts a double value to a string with idigits
    ! decimal places in scientific notation.
!</description>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of 32 characters
!</result>

!<input>

    ! value to be converted
    real(DP), intent(in) :: dvalue

    !number of decimals
    integer, intent(in)  :: idigits
!</input>
!</function>

    soutput = adjustl(sys_sdE(dvalue, idigits))
  end function sys_sdEL


!************************************************************************


!<function>

  function sys_siL(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts an integer value to a string of length idigits,
    ! filled up with white spaces.
!</description>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
!</result>

!<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
!</input>
!</function>

    soutput = adjustl(sys_si(ivalue, idigits))
  end function sys_siL


!************************************************************************


!<function>

  function sys_si0L(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts an integer value to a string of length idigits,
    ! filled up with zeros.
!</description>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
!</result>

!<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
!</input>
!</function>

    soutput = adjustl(sys_si0(ivalue, idigits))
  end function sys_si0L


!************************************************************************


!<function>

  function sys_sliL(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts a long integer value to a string of length idigits.
!</description>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
!</result>

!<input>

    ! value to be converted
    integer(I64), intent(in) :: ivalue

    !number of decimals
    integer, intent(in)      :: idigits
!</input>
!</function>

    soutput = adjustl(sys_sli(ivalue, idigits))
  end function sys_sliL


!************************************************************************


!<function>

  function sys_sli0L(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts a long integer value to a string of length idigits.
!</description>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
!</result>

!<input>

    ! value to be converted
    integer(I64), intent(in) :: ivalue

    !number of decimals
    integer, intent(in)      :: idigits
!</input>
!</function>

    soutput = adjustl(sys_sli0(ivalue, idigits))
  end function sys_sli0L


!************************************************************************
! Wrapper functions to be downward-compatible
! (documentation is omitted in purpose, it would just inflate this
!  file and we do not want them to be used any more)
!
!  sys_i0[3-5]
!  sys_i[1-4,6,8] sys_i64
!  sys_li12
!  sys_s[3,5,6] sys_s1[4,8] sys_s32 sys_s54 sys_s6[1,3] sys_s84
!  sys_d, sys_r

  ! First: int => string

  character (len=3) function sys_i03(ivalue)
    integer, intent(in) :: ivalue
    sys_i03 = trim(sys_si0L(ivalue, 3))
  end function sys_i03

  character (len=4) function sys_i04(ivalue)
    integer, intent(in) :: ivalue
    sys_i04 = trim(sys_si0L(ivalue, 4))
  end function sys_i04

  character (len=5) function sys_i05(ivalue)
    integer, intent(in) :: ivalue
    sys_i05 = trim(sys_si0L(ivalue, 5))
  end function sys_i05

  character (len=1) function sys_i1(ivalue)
    integer, intent(in) :: ivalue
    sys_i1 = trim(sys_siL(ivalue, 1))
  end function sys_i1

  character (len=2) function sys_i2(ivalue)
    integer, intent(in) :: ivalue
    sys_i2 = trim(sys_siL(ivalue, 2))
  end function sys_i2

  character (len=3) function sys_i3(ivalue)
    integer, intent(in) :: ivalue
    sys_i3 = trim(sys_siL(ivalue, 3))
  end function sys_i3

  character (len=4) function sys_i4(ivalue)
    integer, intent(in) :: ivalue
    sys_i4 = trim(sys_siL(ivalue, 4))
  end function sys_i4

  character (len=6) function sys_i6(ivalue)
    integer, intent(in) :: ivalue
    sys_i6 = trim(sys_siL(ivalue, 6))
  end function sys_i6

  character (len=8) function sys_i8(ivalue)
    integer, intent(in) :: ivalue
    sys_i8 = trim(sys_siL(ivalue, 8))
  end function sys_i8

  character (len=10) function sys_i10(ivalue)
    integer, intent(in) :: ivalue
    sys_i10 = trim(sys_siL(ivalue, 10))
  end function sys_i10

  character (len=12) function sys_i12(ivalue)
    integer, intent(in) :: ivalue
    sys_i12 = trim(sys_siL(ivalue, 12))
  end function sys_i12

  character (len=16) function sys_i16(ivalue)
    integer, intent(in) :: ivalue
    sys_i16 = trim(sys_siL(ivalue, 16))
  end function sys_i16

  character (len=64) function sys_i64(ivalue)
    integer, intent(in) :: ivalue
    sys_i64 = trim(sys_siL(ivalue, 64))
  end function sys_i64

  character (len=12) function sys_li12(ivalue)
    integer(I64) :: ivalue
    sys_li12 = trim(sys_sliL(ivalue, 12))
  end function sys_li12


  ! Now: real => string

  character (len=3) function sys_s3(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s3 = trim(sys_sdL(dvalue, 1))
  end function sys_s3

  character (len=5) function sys_s5(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s5 = trim(sys_sdL(dvalue, 2))
  end function sys_s5

  character (len=6) function sys_s6(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s6 = trim(sys_sdL(dvalue, 2))
  end function sys_s6

  character (len=7) function sys_s14(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s14 = trim(sys_sdL(dvalue, 4))
  end function sys_s14

  character (len=11) function sys_s18(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s18 = trim(sys_sdL(dvalue, 8))
  end function sys_s18

  character (len=7) function sys_s32(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s32 = trim(sys_sdL(dvalue, 2))
  end function sys_s32

  character (len=11) function sys_s54(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s54 = trim(sys_sdL(dvalue, 4))
  end function sys_s54

  character (len=10) function sys_s61(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s61 = trim(sys_sdL(dvalue, 1))
  end function sys_s61

  character (len=10) function sys_s63(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s63 = trim(sys_sdL(dvalue, 3))
  end function sys_s63

  character (len=14) function sys_s84(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s84 = trim(sys_sdL(dvalue, 4))
  end function sys_s84

  character (len=16) function sys_d(dvalue)
    real(DP), intent(in) :: dvalue
    sys_d = trim(sys_sdEL(dvalue, 8))
  end function sys_d

  character (len=16) function sys_r(dvalue)
    real(DP), intent(in) :: dvalue
    sys_r = trim(sys_sdL(dvalue, 12))
  end function sys_r

  ! Now: real => string (in scientific notation)

  character (len=9) function sys_s2E(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s2E = trim(sys_sdEL(dvalue, 2))
  end function sys_s2E

  character (len=11) function sys_s4E(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s4E = trim(sys_sdEL(dvalue, 2))
  end function sys_s4E

  character (len=13) function sys_s6E(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s6E = trim(sys_sdEL(dvalue, 2))
  end function sys_s6E

  character (len=17) function sys_s10E(dvalue)
    real(DP), intent(in) :: dvalue
    sys_s10E = trim(sys_sdEL(dvalue, 2))
  end function sys_s10E

end module fsystem
