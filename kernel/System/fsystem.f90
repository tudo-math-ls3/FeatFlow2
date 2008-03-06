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
!#  1.) system_init = system_init_simple /
!#                    system_init_ext
!#      -> Initialise system-wide settings
!#
!#  2.) sys_version
!#      -> Get kernal version number
!#
!#  3.) sys_toupper = sys_toupper_replace /
!#                    sys_toupper_copy
!#      -> Convert a string to uppercase
!#
!#  4.) sys_tolower = sys_tolower_replace /
!#                    sys_tolower_copy
!#      -> Convert a string to lowercase
!#
!#  5.) sys_upcase
!#      -> Convert a string to uppercase, function version
!#
!#  6.) sys_lowcase
!#      -> Convert a string to lowercase, function version
!#
!#  7.) sys_charreplace
!#      -> Replaces characters in a string
!#
!#  8.) sys_throwFPE
!#      -> Throw a floating point exception
!#
!#  9.) sys_getFreeUnit
!#      -> Determine a free file handle for use in an OPEN() command
!#
!# 10.) sys_halt
!#      -> Halts the application. Replacement for CALL sys_halt() in F90.
!#         Can be configured how to halt.
!#         E.g. if the global variable sys_haltmode is set to SYS_HALT_THROFPE,
!#         this routine will stop the program by a floating point exception,
!#         which prints the stack trace to the terminal on some compilers.
!# 
!# 11.) sys_permute
!#      -> Compute a random permutation of a given sequence
!#
!#  ... (documentation incomplete)
!# </purpose>
!##############################################################################

MODULE fsystem

  IMPLICIT NONE

!<constants>

  !<constantblock description="constants for logical values">

  ! logical value 'true'
  INTEGER, PARAMETER :: YES = 0

  ! logical value 'false'
  INTEGER, PARAMETER :: NO = 1

  !</constantblock>
  
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
  INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,37)

  !<constantblock description="Kind values for integers">

  ! kind value for 32Bit integer
  INTEGER, PARAMETER :: I32 = SELECTED_INT_KIND(8)

  ! kind value for 64Bit integer
  INTEGER, PARAMETER :: I64 = SELECTED_INT_KIND(10)

  !</constantblock>

!<constantblock description="system flags">

  ! constant for a system beep
  CHARACTER(len=1), parameter :: BEEP = ACHAR(7)

  ! constant for breaking line in a string
  CHARACTER(len=1), PARAMETER :: NEWLINE = ACHAR(10)

  ! standard length for strings in FEAT
  INTEGER, PARAMETER :: SYS_STRLEN = 256

  ! standard length for name strings in FEAT
  INTEGER, PARAMETER :: SYS_NAMELEN = 32

  ! minimal difference to unity for real values
  REAL(DP)           :: SYS_EPSREAL = EPSILON(1.0_DP)

  ! minimal positive values for real variables
  REAL(DP)           :: SYS_MINREAL = TINY(1.0_DP)

  ! maximal values for real variables
  REAL(DP)           :: SYS_MAXREAL = HUGE(1.0_DP)

  ! maximal values for integer variables
  INTEGER            :: SYS_MAXINT = HUGE(1)

  ! mathematical constant Pi
  REAL(DP)           :: SYS_PI
  
  ! internal constant for infinity
  REAL(DP), PARAMETER :: SYS_INFINITY = HUGE(1.0_DP)

  ! increment value = 1
  INTEGER, PARAMETER :: INCX = 1

  ! flag for appending data to a file
  INTEGER, PARAMETER :: SYS_APPEND = 0

  ! flag for replacing a file
  INTEGER, PARAMETER :: SYS_REPLACE = 1

!</constantblock>
  
!<constantblock description="system signals">

  INTEGER, PARAMETER :: SIGILL  =  4
  INTEGER, PARAMETER :: SIGTRAP =  5
  INTEGER, PARAMETER :: SIGABRT =  6
  INTEGER, PARAMETER :: SIGEMT  =  7
  INTEGER, PARAMETER :: SIGFPE  =  8
  INTEGER, PARAMETER :: SIGBUS  = 10
  INTEGER, PARAMETER :: SIGSEGV = 11
  
!</constantblock>

!<constantblock description="Constants for the sys_haltmode variable">

  ! Halts the program by the CALL sys_halt() command
  INTEGER, PARAMETER :: SYS_HALT_STOP     = 0

  ! Halts the program by sys_throwFPE. On some compilers, this helps with
  ! debugging as the compiler will print a stack trace to the terminal
  ! that allows tracing back where an error came from.
  INTEGER, PARAMETER :: SYS_HALT_THROWFPE = 1
  
!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Emulation of an array of pointers to double precision vectors
  TYPE t_realPointer
    REAL(DP), DIMENSION(:), POINTER :: ptr
  END TYPE

!</typeblock>

!<typeblock>

  ! Emulation of an array of pointers to 2D double precision vectors
  TYPE t_realw2DPointer
    REAL(DP), DIMENSION(:,:), POINTER :: ptr
  END TYPE

!</typeblock>

!<typeblock>

  ! Emulation of an array of pointers to single precision vectors
  TYPE t_singlePointer
    REAL(SP), DIMENSION(:), POINTER :: ptr
  END TYPE

!</typeblock>
  
!<typeblock>

  ! Emulation of an array of pointers to 2D single precision vectors
  TYPE t_single2DPointer
    REAL(SP), DIMENSION(:), POINTER :: ptr
  END TYPE

!</typeblock>

!<typeblock>

  ! Emulation of an array of pointers to integer vectors
  TYPE t_intPointer
    INTEGER(I32), DIMENSION(:), POINTER :: ptr
  END TYPE
  
!</typeblock>
  
!<typeblock>

  ! Emulation of an array of pointers to 2D integer vectors
  TYPE t_int2DPointer
    INTEGER(I32), DIMENSION(:,:), POINTER :: ptr
  END TYPE
  
!</typeblock>
  
!<typeblock>

  ! Global project settings and system information
  TYPE t_sysconfig

    ! project id
    CHARACTER(LEN=SYS_STRLEN) :: sprojectID    = ''

    ! project directory
    CHARACTER(LEN=SYS_STRLEN) :: sprojectDir   = ''

    ! starting time of this project
    INTEGER                   :: iprojectStart = 0

    ! starting time of this project (long time runs)
    INTEGER, DIMENSION(8)     :: iprojectStartLong = 0
    
  END TYPE

!</typeblock>

!<typeblock>

  ! Global time measurement
  TYPE t_clock

    ! The name of the clock
    CHARACTER(LEN=SYS_NAMELEN) :: sname = ''

    ! Starting time
    INTEGER  :: istart       = 0

    ! Stopping time
    INTEGER  :: istop        = 0

    ! Total time elapsed
    INTEGER :: icount        = 0

    ! Total overflow of system clock
    INTEGER :: ioverflow     = 0

    ! Activation flag
    LOGICAL  :: bactive      = .FALSE.

    ! Number of encompassing clock
    INTEGER  :: iencompClock = 0

  END TYPE t_clock

!</typeblock>

!</types>

!<globals>
  ! Global time measurement structure
  TYPE(t_clock), DIMENSION(:), ALLOCATABLE, TARGET, SAVE :: rclock

  ! maximal measurable time span in seconds (system-dependend)
  REAL(DP) :: sys_dtimeMax

  ! global system configuration
  TYPE (t_sysconfig), TARGET, SAVE :: sys_sysconfig
  
  ! Halt mode. This variable defines the way, sys_halt halts the program.
  ! One of the SYS_HALT_xxxx constants.
  INTEGER :: sys_haltmode = SYS_HALT_STOP
!</globals>

  INTERFACE system_init
    MODULE PROCEDURE system_init_simple
    MODULE PROCEDURE system_init_ext
  END INTERFACE

  INTERFACE sys_toupper
    MODULE PROCEDURE sys_toupper_replace
    MODULE PROCEDURE sys_toupper_copy
  END INTERFACE

  INTERFACE sys_tolower
    MODULE PROCEDURE sys_tolower_replace
    MODULE PROCEDURE sys_tolower_copy
  END INTERFACE

CONTAINS

!************************************************************************

!<subroutine>

  SUBROUTINE sys_initClock(nclock)

!<description>

    ! This routine initializes the time measurement.

!</description>

!<input>

    ! Number of clocks.
    INTEGER, INTENT(IN) :: nclock

!</input>

!</subroutine>

    IF (ALLOCATED(rclock)) THEN
      PRINT *, "sys_initClock: Time measurement is already initialised!"
      CALL sys_halt()
    END IF
    
    ALLOCATE(rclock(nclock))
    
  END SUBROUTINE sys_initClock

!************************************************************************

!<subroutine>

  SUBROUTINE sys_doneClock

!<description>

    ! This routine releases the time measurement.

!</description>

!</subroutine>

    IF (.NOT.ALLOCATED(rclock)) THEN
      PRINT *, "sys_doneClock: Time measurement is not initialised!"
      CALL sys_halt()
    END IF
    
    DEALLOCATE(rclock)
    
  END SUBROUTINE sys_doneClock


!************************************************************************

!<subroutine>

  SUBROUTINE sys_setClock(iclock, sname)

!<description>

  ! This routine sets a clock for time measurement.

!</description>

!<input>

    ! Number of the clock
    INTEGER, INTENT(IN)          :: iclock

    ! Name of the clock
    CHARACTER(LEN=*), INTENT(IN) :: sname

!</input>

!</subroutine>

    IF (.NOT.ALLOCATED(rclock)) THEN
      PRINT *, "sys_setClock: Time measurement is not initialised!"
      CALL sys_halt()
    END IF
    
    IF (iclock .GT. SIZE(rclock)) THEN
      PRINT *, "sys_setClock: Clock number exceeds maximum number of clocks!"
      CALL sys_halt()
    END IF
    
    rclock(iclock)%sname        = sname
    rclock(iclock)%istart       = 0
    rclock(iclock)%istop        = 0
    rclock(iclock)%icount       = 0
    rclock(iclock)%ioverflow    = 0
    rclock(iclock)%bactive      = .FALSE.
    rclock(iclock)%iencompClock = 0

  END SUBROUTINE sys_setClock

!************************************************************************

!<subroutine>

  SUBROUTINE sys_setEncompClock(iclock, iencompClock)

!<description>

  ! This routine sets an encompassing clock for a given clock.

!</description>

!<input>

    ! Number of the clock
    INTEGER, INTENT(IN) :: iclock

    ! Number of the encompassing clock
    INTEGER, INTENT(IN) :: iencompClock

!</input>
    
!</subroutine>

    IF (.NOT.ALLOCATED(rclock)) THEN
      PRINT *, "sys_setEncompClock: Time measurement is not initialised!"
      CALL sys_halt()
    END IF
    
    IF (MAX(iclock,iencompClock) .GT. SIZE(rclock)) THEN
      PRINT *, "sys_setEncompClock: Clock number(s) exceed maximum number of clocks!"
      CALL sys_halt()
    END IF
    
    rclock(iclock)%iencompClock = iencompClock
    
  END SUBROUTINE sys_setEncompClock

!************************************************************************

!<subroutine>

  SUBROUTINE sys_startClock(iclock)

!<description>

    ! This routine starts the given clock.

!</description>

!<input>

    ! Number of the clock
    INTEGER, INTENT(IN) :: iclock

!</input>

!</subroutine>

    IF (.NOT.ALLOCATED(rclock)) THEN
      PRINT *, "sys_startClock: Time measurement is not initialised!"
      CALL sys_halt()
    END IF
    
    IF (iclock .GT. SIZE(rclock)) THEN
      PRINT *, "sys_startClock: Clock number exceeds maximum number of clocks!"
      CALL sys_halt()
    END IF
    
    ! Start time measurement
    rclock(iclock)%bactive = .TRUE.
    CALL system_clock(rclock(iclock)%istart)
  END SUBROUTINE sys_startClock

!************************************************************************

!<subroutine>

  SUBROUTINE sys_stopClock(iclock)

!<description>

    ! This routine stops the given clock.

!</description>

!<input>

    ! Number of the clock
    INTEGER, INTENT(IN) :: iclock

!</input>

!</subroutine>

    ! local variables
    INTEGER :: icount,irate,icountmax

    IF (.NOT.ALLOCATED(rclock)) THEN
      PRINT *, "sys_stopClock: Time measurement is not initialised!"
      CALL sys_halt()
    END IF
    
    IF (iclock .GT. SIZE(rclock)) THEN
      PRINT *, "sys_stopClock: Clock number exceeds maximum number of clocks!"
      CALL sys_halt()
    END IF
    
    IF (.NOT.rclock(iclock)%bactive) THEN
      PRINT *, "sys_stopClock: Clock has not been started, skipping!"
      RETURN
    END IF

    ! Stop time measurement
    rclock(iclock)%bactive = .FALSE.
    CALL system_clock(rclock(iclock)%istop, irate, icountmax)
    
    ! Calculate elapsed time. Note that the stopping time may be smaller
    ! than the starting time. In this case the maximum number of counts
    ! has been reached and must be considered.
    IF (rclock(iclock)%istart  .GT. rclock(iclock)%istop) THEN
      icount = icountmax-rclock(iclock)%istart
      icount = icount+rclock(iclock)%istop
      IF (rclock(iclock)%icount .GT. icountmax-icount) THEN
        rclock(iclock)%icount    = rclock(iclock)%icount-icountmax
        rclock(iclock)%icount    = rclock(iclock)%icount+icount
        rclock(iclock)%ioverflow = rclock(iclock)%ioverflow+1
      ELSE
        rclock(iclock)%icount =  rclock(iclock)%icount+icount
      END IF
    ELSE
      icount = rclock(iclock)%istop-rclock(iclock)%istart
      IF (rclock(iclock)%icount .GT. icountmax-icount) THEN
        rclock(iclock)%icount    = rclock(iclock)%icount-icountmax
        rclock(iclock)%icount    = rclock(iclock)%icount+icount
        rclock(iclock)%ioverflow = rclock(iclock)%ioverflow+1
      ELSE
        rclock(iclock)%icount = rclock(iclock)%icount+icount
      END IF
    END IF
  END SUBROUTINE sys_stopClock

!************************************************************************

!<subroutine>

  SUBROUTINE sys_stopClockAll

!<description>

    ! This routine stops all clocks.

!</description>

!</subroutine>

    ! local variables
    INTEGER :: iclock

    IF (.NOT.ALLOCATED(rclock)) THEN
      PRINT *, "sys_stopClockAll: Time measurement is not initialised!"
      CALL sys_halt()
    END IF

    DO iclock = 1, SIZE(rclock)
      CALL sys_stopClock(iclock)
    END DO

  END SUBROUTINE sys_stopClockAll

!************************************************************************

!<subroutine>

  SUBROUTINE sys_infoClock

!<description>

    ! This routine prints information about the clock

!</description>

!</subroutine>

    ! local variables
    TYPE(t_clock), DIMENSION(:), POINTER :: rclockTmp
    REAL(DP) :: dtotaltime,dtime,dtimeself
    INTEGER  :: iclock,icount,irate,icountmax
    INTEGER  :: ndays,nsecs,nsec1,nsec2
    INTEGER, DIMENSION(8) :: istopLong

    WRITE(*,FMT='(A)') 'Time measurement:'
    WRITE(*,FMT='(A)') '-----------------'

    ! Compute total cpu time
    CALL system_clock(icount,irate,icountmax)
    CALL date_and_time(values=istopLong)

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

    IF (nsec2 .LT. nsec1) THEN
      nsecs = 86400-(nsec2-nsec1)
    ELSE
      nsecs = nsec2-nsec1
    END IF

    ! Compute total time from system clock
    IF (icount .LT. sys_sysconfig%iprojectStart) THEN
      dtime = REAL(icountmax+icount-sys_sysconfig%iprojectStart, DP)/REAL(irate, DP)
    ELSE
      dtime = REAL(icount-sys_sysconfig%iprojectStart, DP)/REAL(irate, DP)
    END IF

    ! Ok, now we have two times. Check which measurement can be used
    dtotaltime = REAL(icountmax, DP)/REAL(irate, DP)

    IF (ndays .LE. INT(dtotaltime/86400._DP)) THEN
      dtotaltime = dtime
      WRITE(*,FMT='(A)') 'Simulation ran '//TRIM(sys_sdEL(dtime, 2))//' seconds'
    ELSE
      dtotaltime = REAL(ndays*86400+nsecs, DP)
      WRITE(*,FMT='(A)') 'Simulation ran '//TRIM(sys_siL(ndays, 15))//' days and '//&
                                            TRIM(sys_siL(nsecs, 15))//' seconds'
    END IF

    IF (.NOT.ALLOCATED(rclock)) THEN
      PRINT *, "sys_infoClock: Time measurement is not initialised!"
      CALL sys_halt()
    END IF

    ! Make a copy of the system clock
    ALLOCATE(rclockTmp(SIZE(rclock)))
    rclockTmp = rclock

    ! Subtract the individual times from the encompassing clocks
    DO iclock = 1, SIZE(rclockTmp)
      IF (rclockTmp(iclock)%iencompClock .NE. 0) THEN
        rclockTmp(rclockTmp(iclock)%iencompClock)%icount =&
            rclockTmp(rclockTmp(iclock)%iencompClock)%icount-rclockTmp(iclock)%icount
        IF (rclockTmp(rclockTmp(iclock)%iencompClock)%icount .LT. 0) THEN
          rclockTmp(rclockTmp(iclock)%iencompClock)%icount =&
              rclockTmp(rclockTmp(iclock)%iencompClock)%icount+icountmax
          rclockTmp(rclockTmp(iclock)%iencompClock)%ioverflow =&
              rclockTmp(rclockTmp(iclock)%iencompClock)%ioverflow-1  
        END IF
      END IF
    END DO
    
    ! Check consistency
    DO iclock = 1, SIZE(rclockTmp)
      IF (rclockTmp(iclock)%icount .LT. 0 .OR.&
          rclockTmp(iclock)%ioverflow .LT. 0) THEN
        PRINT *, "sys_infoClock: Time measurement is incorrect!"
        CALL sys_halt()
      END IF
    END DO

    WRITE(*,*)
    WRITE(*,FMT='(A,T35,A,T60,A)') 'Clock name', 'elapsed time (%)', 'self time (%)'
    WRITE(*,FMT='(85("="))')

    ! Print out all clocks
    DO iclock = 1, SIZE(rclock)
      
      ! Compute time in seconds
      dtime = rclock(iclock)%ioverflow*REAL(icountmax, DP)/REAL(irate, DP)+&
              REAL(rclock(iclock)%icount, DP)/REAL(irate, DP)
      dtimeself = rclockTmp(iclock)%ioverflow*REAL(icountmax, DP)/REAL(irate, DP)+&
                  REAL(rclockTmp(iclock)%icount, DP)/REAL(irate, DP)

      WRITE(*,FMT='(A,T35,A,T60,A)')&
          TRIM(ADJUSTL(rclock(iclock)%sname)),&
          TRIM(sys_sdEL(dtime, 2))//' ('//&
          TRIM(sys_sdL(100._DP/dtotaltime*dtime, 2))//')',&
          TRIM(sys_sdEL(dtimeself, 2))//' ('//&
          TRIM(sys_sdL(100._DP/dtotaltime*dtimeself,2))//')'
    END DO

    ! Free unused memory
    DEALLOCATE(rclockTmp)

  CONTAINS

    FUNCTION calender_to_julian(year, month, day) RESULT(ivalue)
      INTEGER, INTENT(IN) :: year
      INTEGER, INTENT(IN) :: month
      INTEGER, INTENT(IN) :: day

      INTEGER             :: ivalue

      ivalue = day-32075+&
               1461*(year+4800+(month-14)/12)/4+&
               367*(month-2-((month-14)/12)*12)/12-&
               3*((year+4900+(month-14)/12)/100)/4
    END FUNCTION calender_to_julian
  END SUBROUTINE sys_infoClock

!************************************************************************

!<subroutine>

  SUBROUTINE sys_permute(k,Idata)

!<description>
    ! This routine computes the permutation of the initial set Idata
    ! that corresponds to the factoriodic number k.
!</description>

!<input>
    ! factoriodic number
    INTEGER(I32), INTENT(IN) :: k
!</input>

!<inputoutput>
    ! initial and permuted set
    INTEGER(I32), DIMENSION(:), INTENT(INOUT) :: Idata
!</inputoutput>
!</subroutine>
  
    ! local variables
    INTEGER(I32) :: factorial,i,j,iswap
    
    DO j=2,SIZE(Idata)

      IF (j-1 .LT. k) THEN
        factorial = 0._I64
      ELSEIF(j-1 .EQ. k) THEN
        factorial = 1._I64
      ELSE
        factorial=k+1
        DO i=k+2,j-1
          factorial=factorial*i
        END DO
      END IF
      
      i = j-MOD(factorial,j)
      iswap    = Idata(j)
      Idata(j) = Idata(i)
      IData(i) = iswap
    END DO
  END SUBROUTINE sys_permute

!************************************************************************

!<subroutine>

  SUBROUTINE sys_halt()
  
!<description>
  ! This routine halts the application like the CALL sys_halt() command in
  ! Fortran 90. The routine can be configured how to halt the application.
  ! For this purpose, the main program can set the global variable
  ! sys_haltmode to one of the SYS_HALT_xxxx constants.
!</description>
    
!</subroutine>

    SELECT CASE (sys_haltmode)
    CASE (SYS_HALT_STOP)
      STOP
    CASE (SYS_HALT_THROWFPE)
      CALL sys_throwFPE()
    END SELECT

  END SUBROUTINE sys_halt

!************************************************************************

!<subroutine>

  PURE SUBROUTINE sys_throwFPE()
  
    !<description>
    !This routine throws a floating point exception for debugging purposes
    !to prevent the debugger to exit the program.
    !</description>
    
!</subroutine>

    INTEGER :: i1,i2

    i1=1
    i2=0
    i1=i1/i2
    
  END SUBROUTINE sys_throwFPE

!************************************************************************

!<subroutine>

  SUBROUTINE system_init_simple()

!<description>
  ! This subroutine initialises internal data structures with standard
  ! values.
!</description>

!</subroutine>

    CALL system_init_ext("","")

  END SUBROUTINE 

!************************************************************************

!<subroutine>

  SUBROUTINE system_init_ext(sprojectID,sprojectDir)

!<description>
  ! Extended initialisation.
  ! This subroutine initialises internal data structures.
  ! The value of the project name and directory are set according to the
  ! parameters.
!</description>

!<input>
  ! An ID string of the project.
  CHARACTER(LEN=*), INTENT(IN) :: sprojectID

  ! The directory of the project. "" means 'current directory'.
  CHARACTER(LEN=*), INTENT(IN) :: sprojectDir
!</input>

!</subroutine>

    INTEGER :: icount ! current system time
    INTEGER :: irate  ! approx. number of system clock ticks per second
    INTEGER :: icmax  ! largest possible value of icount

    ! system_clock is not a FEAT, but a basic FORTRAN 90 routine
    CALL system_clock(icount,irate,icmax)

    ! use data_and_time to measure long time runs
    CALL date_and_time(values=sys_sysconfig%iprojectStartLong)

    ! maximal measurable time span in seconds (system-dependend)
    sys_dtimeMax = REAL(icmax,DP)/REAL(irate,DP)

    ! Initialise the global sysconfig structure
    sys_sysconfig%sprojectID    = sprojectID
    sys_sysconfig%sprojectDir   = sprojectDir
    sys_sysconfig%iprojectStart = icount

    ! Set value of Pi = 3.14..
    SYS_PI=ASIN(1.0_DP)*2.0_DP

  END SUBROUTINE

!************************************************************************************


!<subroutine>
  SUBROUTINE sys_version(ifeastVersionHigh, ifeastVersionMiddle, ifeastVersionLow, &
                         sreldate)

!<description>
  ! This subroutine returns the library version information.
!</description>

!<output>

  ! high version number
  INTEGER :: ifeastVersionHigh

  ! middle version number
  INTEGER :: ifeastVersionMiddle

  ! low version number
  INTEGER :: ifeastVersionLow

  ! release date
  CHARACTER(LEN=*) :: sreldate

!</output>

!</subroutine>

    ifeastVersionHigh=0
    ifeastVersionMiddle=0
    ifeastVersionLow=1

    sreldate="01.01.2007 RC0"

  END SUBROUTINE sys_version

!************************************************************************

!<subroutine>
  
  SUBROUTINE sys_toupper_replace (str)

!<description>
  ! Convert a string to upper case. 
  ! The given string is replaced by its uppercase version.
!</description>

!<inputoutput>
  
  ! The string that is to make uppercase
  CHARACTER(LEN=*), INTENT(INOUT) :: str
  
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  
  INTEGER, PARAMETER :: up2low = IACHAR("a") - IACHAR("A")
  INTEGER :: i
  CHARACTER    :: c
      
  DO i=1,LEN(str)
    c = str(i:i)
    IF ((c .GE. "a") .AND. (c .LE. "z")) THEN
      str(i:i) = ACHAR (IACHAR(c) - up2low)
    END IF
  END DO
  END SUBROUTINE

!************************************************************************

!<subroutine>
  
  SUBROUTINE sys_toupper_copy (str,strUpper) 

!<description>
  ! Convert a string to upper case.
!</description>

!<input>
  
  ! The string that is to make uppercase
  CHARACTER(LEN=*), INTENT(IN) :: str

!</input>

!<output>

  ! Uppercase version of the given string
  CHARACTER(LEN=*), INTENT(OUT) :: strUpper
  
!</output>
  
!</subroutine>
  
  ! local variables
  
  INTEGER, PARAMETER :: up2low = IACHAR("a") - IACHAR("A")
  INTEGER :: i
  CHARACTER    :: c

  IF (LEN(str) > LEN(strUpper)) THEN
    PRINT *, "sys_toupper_copy: target string is too short"
    CALL sys_halt()
  END IF
  
  ! Initialize string
  strUpper = ''
  
  DO i=1,LEN(str)
    c = str(i:i)
    IF ((c .GE. "a") .AND. (c .LE. "z")) THEN
      strUpper(i:i) = ACHAR (IACHAR(c) - up2low)
    ELSE
      strUpper(i:i) = c
    END IF
  END DO
  END SUBROUTINE

!************************************************************************

!<subroutine>
  
  SUBROUTINE sys_tolower_replace (str) 

!<description>
  ! Convert a string to lower case.
  ! The given string is replaced by its lowercase version.
!</description>

!<inputoutput>
  
  ! The string that is to make lowercase
  CHARACTER(LEN=*), INTENT(INOUT) :: str

!</inputoutput>
  
!</subroutine>
  
  ! local variables
  
  INTEGER, PARAMETER :: up2low = IACHAR("a") - IACHAR("A")
  INTEGER :: i
  CHARACTER    :: c
  
  DO i=1,LEN(str)
    c = str(i:i)
    IF ((c .GE. "A") .AND. (c .LE. "Z")) THEN
      str(i:i) = ACHAR (IACHAR(c) + up2low)
    END IF
  END DO
  END SUBROUTINE
  
!************************************************************************

!<subroutine>
  
  SUBROUTINE sys_tolower_copy (str,strLower) 

!<description>
  ! Convert a string to lower case.
!</description>

!<input>
  
  ! The string that is to make lowercase
  CHARACTER(LEN=*), INTENT(IN) :: str

!</input>

!<output>

  ! Lowercase version of the given string
  CHARACTER(LEN=*), INTENT(OUT) :: strLower
  
!</output>
  
!</subroutine>
  
  ! local variables
  
  INTEGER, PARAMETER :: up2low = IACHAR("a") - IACHAR("A")
  INTEGER :: i
  CHARACTER    :: c

  IF (LEN(str) > LEN(strLower)) THEN
    PRINT *, "sys_tolower_copy: target string is too short"
    CALL sys_halt()
  END IF
  
  ! Initialize string
  strLower = ''
  
  DO i=1,LEN(str)
    c = str(i:i)
    IF ((c .GE. "A") .AND. (c .LE. "Z")) THEN
      strLower(i:i) = ACHAR (IACHAR(c) + up2low)
    ELSE
      strLower(i:i) = c
    END IF
  END DO
  END SUBROUTINE

!******************************************************************************

!<function>
  PURE FUNCTION sys_upcase(sinput) RESULT(soutput)
  
  !<description>
  ! This routine converts a given string to its uppercase version.
  !</description>

  !<input>

  !input string
  character(len=*), intent(in) :: sinput
  !</input>

  !<output>

  !output string
  character(len=len(sinput)) :: soutput

  !</output>
!</function>

  !index variable
  INTEGER :: i

  soutput = " "   !initialise string
  DO I = 1,LEN(sinput)
     if(sinput(i:i) .GE. "a" .and. sinput(i:i) .LE. "z") then
        soutput(i:i) = achar(iachar(sinput(i:i)) - 32)
     ELSE
        soutput(i:i) = sinput(i:i)
     END IF
  END DO
  
  END FUNCTION sys_upcase

!******************************************************************************

!<function>
  PURE FUNCTION sys_lowcase(sinput) RESULT(soutput)
  
  !<description>
  ! This routine converts a given string to its uppercase version.
  !</description>

  !<input>

  !input string
  character(len=*), intent(in) :: sinput
  !</input>

  !<output>

  !output string
  character(len=len(sinput)) :: soutput

  !</output>
!</function>

  !index variable
  INTEGER :: i

  soutput = " "   !initialise string
  DO I = 1,LEN(sinput)
     if(sinput(i:i) .GE. "A" .and. sinput(i:i) .LE. "Z") then
        soutput(i:i) = achar(iachar(sinput(i:i)) + 32)
     ELSE
        soutput(i:i) = sinput(i:i)
     END IF
  END DO
  
  END FUNCTION sys_lowcase

!******************************************************************************

!<function>
  PURE FUNCTION sys_charreplace(sinput,scharsource,schardest) RESULT(soutput)
  
!<description>
  ! Replaces all characers scharsource in sinput by schardest.
  ! Case sensitive.
!</description>

!<input>
  ! input string
  CHARACTER(LEN=*), INTENT(IN) :: sinput
  
  ! Character to be searched for.
  CHARACTER, INTENT(IN) :: scharsource
  
  ! Detinatiion character, all scarsource characters in sinput should be
  ! replaced by.
  CHARACTER, INTENT(IN) :: schardest
!</input>

!<output>
  ! output string
  CHARACTER(LEN=LEN(sinput)) :: soutput
!</output>
!</function>

  !index variable
  INTEGER :: i

  soutput = " "   !initialise string
  DO I = 1,LEN(sinput)
     if(sinput(i:i) .eq. scharsource) then
        soutput(i:i) = schardest
     ELSE
        soutput(i:i) = sinput(i:i)
     END IF
  END DO
  
  END FUNCTION 

!************************************************************************

!<function>
  INTEGER function sys_getFreeUnit()

    !<description>
    !This routine tries to find a free unit (for file input/output). If a free unit is
    !found, it is returned, otherwise -1 is returned.
    !</description>

    !<result>
    !number of free unit (-1 if no free unit available)
    !</result>

!</function>

    logical :: bexists, bopened!flags indicating errors
    INTEGER :: itry !free unit candidate

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
    INTEGER :: iunit

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

!<function>
  subroutine sys_flush(iunit)

    !<description>
    ! This routine flushes the buffers associated with an open output unit.
    ! This normally happens when the file is closed or the program ends, 
    ! but this routine ensures the buffers are flushed before any other 
    ! processing occurs.
    !</description>

    !<input>

    !unit connected to the file to write to
    INTEGER :: iunit

    !</input>
!</function>

!#ifdef HAS_FLUSH
!    call flush(iunit)
!#endif

  end subroutine sys_flush


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
    logical, INTENT(IN) :: lvalue
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

END MODULE
