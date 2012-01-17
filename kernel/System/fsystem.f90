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
!#  2.) sys_halt
!#      -> Halts the application. Replacement for STOP commands in F90.
!#         Can be configured how to halt.
!#         E.g. if the global variable sys_haltmode is set to SYS_HALT_THROWFPE,
!#         this routine will stop the program by a floating point exception,
!#         which prints the stack trace to the terminal on some compilers.
!#
!#  3.) sys_version
!#      -> Get kernel version number
!#
!#  4.) sys_throwFPE
!#      -> Throw a floating point exception
!#
!#  5.) sys_permute
!#      -> Compute a random permutation of a given sequence
!#
!#  6.) sys_toupper = sys_toupper_replace /
!#                    sys_toupper_copy
!#      -> Convert a string to uppercase
!#
!#  7.) sys_tolower = sys_tolower_replace /
!#                    sys_tolower_copy
!#      -> Convert a string to lowercase
!#
!#  8.) sys_upcase
!#      -> Convert a string to uppercase, function version
!#
!#  9.) sys_lowcase
!#      -> Convert a string to lowercase, function version
!#
!# 10.) sys_charreplace
!#      -> Replaces characters in a string
!#
!# 11.) sys_getFreeUnit
!#      -> Determine a free file handle for use in an OPEN() command
!#
!# 12.) sys_fileExists
!#      -> Check if file with a given name exists
!#
!# 13.) sys_flush
!#      -> Flush file (if available)
!#
!# 14.) sys_str2Double
!#      -> Convert string to double value
!#
!# 15.) sys_str2Single
!#      -> Convert string to single value
!#
!# 16.) sys_d   ,sys_sd  ,sys_sdP  ,sys_sdE ,sys_sdEP ,sys_sdL ,sys_sdEL,
!#      sys_r   ,sys_s3  ,sys_s5   ,sys_s6  ,sys_s14  ,sys_s18 ,sys_s32,
!#      sys_s54 ,sys_s61 ,sys_s63  ,sys_s84 ,
!#      sys_s2E ,sys_s4E ,sys_s6E  ,sys_s10E
!#      -> String routines to convert double precision numbers to strings
!#
!# 17.) sys_si  ,sys_si0 ,sys_sli ,sys_sli0 ,sys_siL ,sys_si0L ,
!#      sys_i03 ,sys_i04 ,sys_i05 ,sys_i1   ,sys_i2  ,sys_i3   ,sys_i4,
!#      sys_i6  ,sys_i8  ,sys_i10 ,sys_i12  ,sys_i16 ,sys_i64
!#      -> String routines to convert integer numbers to strings
!#
!# 18.) sys_sliL, sys_sli0L, sys_li12
!#      -> String routines to convert long integer numbers to strings
!#
!# 19.) sys_sl
!#      -> String routines to convert a logical to a strings
!#
!# 20.) sys_smem, sys_smemL
!#      -> String routines to convert long integers representing memory usage
!#         to strings.
!#
!# 21.) sys_getenv_string
!#      -> Retrieves an environment variable from the system
!#
!# 22.) sys_silsb, sys_simsb
!#      -> String routine to convert integer numbers to string in bis
!#         representation (LSB/MSB)
!#
!# 23.) sys_adjustl, sys_adjustr
!#      -> Extended version of ADJUSTL/ADJUSTR that allow to specify
!#         the length of the resulting string as parameter
!#
!# 24.) sys_dequote
!#      -> De-quote a string; remove any quotation marks
!#
!# 25.) sys_stringToCharArray
!#      -> Converts a string into a character array
!#
!# 26.) sys_charArrayToString
!#      -> Converts a character array into a string
!#
!# </purpose>
!##############################################################################

module fsystem

!$use omp_lib

  implicit none

!<constants>

!<constantblock description="constants for logical values">

  ! logical value 'true'
  integer, parameter :: YES = 0

  ! logical value 'false'
  integer, parameter :: NO = 1

!</constantblock>
 
!<constantblock description="Kind values for floats">

  ! kind value for 32 bit float (single precision)
  integer, parameter :: SP = selected_real_kind(6,37)

  ! kind value for 64 bit float (double precision)
  integer, parameter :: DP = selected_real_kind(15,307)

#ifdef ENABLE_QUADPREC
  ! kind value for 80/128 bit float (quad precision)
  integer, parameter :: QP = selected_real_kind(18,4931)
#else
  ! set QP equal to DP to avoid compiler problems
  integer, parameter :: QP = DP
#endif
  
  ! Note: Depending on the platform and the compiler, QP is either an 80
  ! or an 128 bit float. The g95 and gfortran compilers use 80 floats
  ! for QP, while the ifc compiler uses 128 bit floats.

!</constantblock>

!<constantblock description="Kind values for integers">

  ! kind value for 8 bit integer
  integer, parameter :: I8 = selected_int_kind(2)

  ! kind value for 16 bit integer
  integer, parameter :: I16 = selected_int_kind(4)

  ! kind value for 32 bit integer
  integer, parameter :: I32 = selected_int_kind(8)

  ! kind value for 64 bit integer
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
  real(SP), parameter :: SYS_EPSREAL_SP = epsilon(1.0_SP)
  real(DP), parameter :: SYS_EPSREAL_DP = epsilon(1.0_DP)
  real(QP), parameter :: SYS_EPSREAL_QP = epsilon(1.0_QP)

  ! minimal positive values for real variables
  real(SP), parameter :: SYS_MINREAL_SP = tiny(1.0_SP)
  real(DP), parameter :: SYS_MINREAL_DP = tiny(1.0_DP)
  real(QP), parameter :: SYS_MINREAL_QP = tiny(1.0_QP)

  ! maximal values for real variables
  real(SP), parameter :: SYS_MAXREAL_SP = huge(1.0_SP)
  real(DP), parameter :: SYS_MAXREAL_DP = huge(1.0_DP)
  real(QP), parameter :: SYS_MAXREAL_QP = huge(1.0_QP)

  ! maximal values for integer variables
  integer,      parameter :: SYS_MAXINT = huge(1)
  integer(I8),  parameter :: SYS_MAXI8  = huge(1_I8)
  integer(I16), parameter :: SYS_MAXI16 = huge(1_I16)
  integer(I32), parameter :: SYS_MAXI32 = huge(1_I32)
  integer(I64), parameter :: SYS_MAXI64 = huge(1_I64)

  ! mathematical constant Pi
  real(DP)           :: SYS_PI
  
  ! internal constant for infinity
  real(SP), parameter :: SYS_INFINITY_SP = huge(1.0_SP)
  real(DP), parameter :: SYS_INFINITY_DP = huge(1.0_DP)
  real(QP), parameter :: SYS_INFINITY_QP = huge(1.0_QP)

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
  type t_real2DPointer
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

    ! Starting time of this project (long time runs).
    ! Format: year / month / day / time-difference to UTC /
    !         hours / minutes / seconds / milliseconds.
    integer, dimension(8)     :: iprojectStartLong = 0
    
  end type

!</typeblock>

!</types>

!************************************************************************

!<publicvars>
  ! global system configuration
  type (t_sysconfig), target, save :: sys_sysconfig
  
  ! Halt mode. This variable defines the way, sys_halt halts the program.
  ! One of the SYS_HALT_xxxx constants.
  integer, save :: sys_haltmode = SYS_HALT_STOP

  ! The Fortran system_clock timer, like all integer timers, has a cycle
  ! time of real(max)/real(rate) seconds. After max clock cycles the
  ! clock will start counting again from zero. This is the maximum time
  ! span that can be measured when using system_clock manually.
  !
  ! Note: Timing routines in the statistics module automatically
  ! respect this setting but do not explicitely use this variable.
  real(DP), save :: sys_dtimeMax = 0.0_DP

!</publicvars>

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

  subroutine sys_permute(k,Idata)

!<description>
    ! This routine computes the permutation of the initial set Idata
    ! that corresponds to the factoriodic number k.
!</description>

!<input>
    ! factoriodic number
    integer(I32), intent(in) :: k
!</input>

!<inputoutput>
    ! initial and permuted set
    integer(I32), dimension(:), intent(inout) :: Idata
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
    ! The command line parameters are read and stored into
    ! sys_scommandLineArgs.
!</description>

!<input>
    ! An ID string of the project.
    character(LEN=*), intent(in) :: sprojectID
    
    ! The directory of the project. "" means 'current directory'.
    character(LEN=*), intent(in) :: sprojectDir
!</input>

!</subroutine>

    ! local variables
    integer :: icount ! current system time
    integer :: irate  ! approx. number of system clock ticks per second
    integer :: icmax  ! largest possible value of icount

    ! system_clock is not a FEAT, but a basic FORTRAN 90 routine
    call system_clock(icount,irate,icmax)

    ! compute maximam measurable timespan
    sys_dtimeMax = real(icmax,DP)/real(irate,DP)

    ! use data_and_time to measure long time runs
    call date_and_time(values=sys_sysconfig%iprojectStartLong)

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
    
    sreldate="01.01.2009 RC0"
    
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
    character(LEN=*), intent(inout) :: str
  
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
    character(LEN=*), intent(in) :: str

!</input>

!<output>

    ! Uppercase version of the given string
    character(LEN=*), intent(out) :: strUpper
  
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
    
    ! Initialise string
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
    character(LEN=*), intent(inout) :: str

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
    character(LEN=*), intent(in) :: str

!</input>

!<output>

    ! Lowercase version of the given string
    character(LEN=*), intent(out) :: strLower
  
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
    
    ! Initialise string
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
    character(LEN=*), intent(in) :: sinput
    
    ! Character to be searched for.
    character, intent(in) :: scharsource
    
    ! Detinatiion character, all scarsource characters in sinput should be
    ! replaced by.
    character, intent(in) :: schardest
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

#ifdef HAS_FLUSH
    call flush(iunit)
#endif

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
    character(LEN=*), intent(in) :: svalue

    ! format description to use for conversion
    character(LEN=*), intent(in) :: sformat

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
    character(LEN=*), intent(in) :: svalue

    ! format description to use for conversion
    character(LEN=*), intent(in) :: sformat

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

!<function>
  
  character (len=32) function sys_smem(imem,ndigits) result(sout)

!<description>
  ! This routine converts 64 bit integer representing memory usage to
  ! a string.
!</description>

!<input>
  ! The memory usage in bytes.
  integer(I64), intent(in) :: imem
  
  ! OPTIONAL: The number of digits. Must be 0 <= ndigits <= 3.
  ! If not given, ndigits = 2 is used.
  integer, optional, intent(in) :: ndigits
!</input>

!</function>

  integer :: k,nld,ntd
  integer(I64) :: ii,jj,itds
  character(len=32) :: sformat
  character(len=8) :: spost
  character(len=2) :: snld, sntd
  
    ! Get the number of trailing digits
    ntd = 2
    if(present(ndigits)) then
      if(ndigits .lt. 0) then
        ntd = 0
      else if(ndigits .gt. 3) then
        ntd = 3
      else
        ntd = ndigits
      end if
    end if
    
    ! Calculate trailing digits scale
    itds = (10_I64)**ntd
  
    ! Find out in which range the memory usage is:
    ii = abs(imem)
    jj = 0
    do k = 0, 6
      if(ii .lt. 1024_I64) exit
      jj = (itds * mod(ii, 1024_I64)) / 1024_I64
      ii = ii / 1024_I64
    end do
  
    ! What do we have here?
    select case(k)
    case (0)
      spost = ' Bytes'
    case (1)
      spost = ' KB'   ! "Kilobytes"
    case (2)
      spost = ' MB'   ! "Megabytes"
    case (3)
      spost = ' GB'   ! "Gigabytes"
    case (4)
      spost = ' TB'   ! "Terabytes"
    case (5)
      spost = ' PB'   ! "Petabytes"
    case (6)
      spost = ' EB'   ! "Exabytes"
    end select
    
    ! "Count" the number of leading digits
    if(ii .lt. 10_I64) then
      nld = 1
    else if(ii .lt. 100_I64) then
      nld = 2
    else if(ii .lt. 1000_I64) then
      nld = 3
    else
      nld = 4
    end if
    
    ! If the memory usage was negative (nice idea), then the number
    ! of leading digits has to be increased by 1 for the sign.
    if(imem .lt. 0_I64) then
      nld = nld + 1
      ii = -ii
    end if
    
    ! Prepare snld and sntd
    write(snld,'(i1)') nld
    write(sntd,'(i1)') ntd
    
    ! Now what format are we going to print?
    if((k .eq. 0) .or. (ntd .eq. 0)) then

      ! Print something like "xxx KB"
      sformat = '(i' // trim(snld) // ',"' // trim(spost) // '")'
      write(sout, sformat) ii
      
    else

      ! Print something like "xxx.yy KB"
      sformat = '(i' // trim(snld) // ',".",i' // trim(sntd) // '.' &
                // trim(sntd) // ',"' // trim(spost) // '")'
      write(sout, sformat) ii, jj
      
    end if
    
  end function sys_smem

!************************************************************************

!<function>
  
  character (len=32) function sys_smemL(imem,ndigits) result(sout)

!<description>
  ! This routine converts 64 bit integer representing memory usage to
  ! a string.
!</description>

!<input>
  ! The memory usage in bytes.
  integer(I64), intent(in) :: imem
  
  ! OPTIONAL: The number of digits. Must be 0 <= ndigits <= 3.
  ! If not given, ndigits = 2 is used.
  integer, optional, intent(in) :: ndigits
!</input>

!</function>
  
    if(present(ndigits)) then
      sout = adjustl(sys_smem(imem,ndigits))
    else
      sout = adjustl(sys_smem(imem))
    end if

  end function

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
    logical, intent(in) :: lvalue
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

    if (ipositions .gt. 32) then
      write(6, *) "*** WARNING! Too long string requested in sys_sdP! ***"
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

    if (ipositions .gt. 24) then
      write(6, *) "*** WARNING! Too long string requested in sys_sdEP! ***"
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

!<input>
    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
!</input>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
!</result>

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

!<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
!</input>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
!</result>
!</function>

    soutput = adjustl(sys_si0(ivalue, idigits))
  end function sys_si0L


!************************************************************************


!<function>

  function sys_sliL(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts a long integer value to a string of length idigits.
!</description>

!<input>
    ! value to be converted
    integer(I64), intent(in) :: ivalue

    !number of decimals
    integer, intent(in)      :: idigits
!</input>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
!</result>
!</function>

    soutput = adjustl(sys_sli(ivalue, idigits))
  end function sys_sliL


!************************************************************************


!<function>

  function sys_sli0L(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts a long integer value to a string of length idigits.
!</description>

!<input>
    ! value to be converted
    integer(I64), intent(in) :: ivalue

    !number of decimals
    integer, intent(in)      :: idigits
!</input>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
!</result>
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

!************************************************************************


!<function>

  function sys_adjustr(sstring,nchars) result(soutput)

!<description>
    ! Extended ADJUSTR.
    ! Right-adjusts a string onto a length on nchars characters
!</description>

!<input>
    ! value to be converted
    character(len=*), intent(in) :: sstring

    ! number of characters
    integer, intent(in)      :: nchars
!</input>

!<result>
    ! String representation of the value (right-aligned),
    ! fixed length of nchars characters
    character (len=nchars) :: soutput
!</result>
!</function>

    write(soutput,"(A"//trim(sys_siL(nchars,10))//")") trim(adjustl(sstring))
    
  end function

!************************************************************************


!<function>

  function sys_adjustl(sstring,nchars) result(soutput)

!<description>
    ! Extended ADJUSTR.
    ! Left-adjusts a string onto a length on nchars characters
!</description>

!<input>
    ! value to be converted
    character(len=*), intent(in) :: sstring

    ! number of characters
    integer, intent(in)      :: nchars
!</input>

!<result>
    ! String representation of the value (right-aligned),
    ! fixed length of nchars characters
    character (len=nchars) :: soutput
!</result>
!</function>

    soutput = trim(adjustl(sstring))
    
  end function

  ! ***************************************************************************

!<function>
  logical function sys_getenv_string(svar, sresult)

  !<description>
    ! This functions returns the string value of a given enviroment variable. The routine
    ! returns .TRUE., if the variable exists, otherwise .FALSE. .
  !</description>

  !<input>
    ! name of the enviroment variable
    character(len=*), intent(in) :: svar
  !</input>

  !<output>
    ! value of the enviroment variable
    character(len=SYS_STRLEN), intent(out) :: sresult
  !</output>

  !<result>
    ! exit status
  !</result>

  !<errors>
    ! none
  !</errors>
!</function>

    character(len=SYS_STRLEN) :: svalueInEnv

    integer :: nstatus

#if defined USE_COMPILER_NEC || defined USE_COMPILER_PGI_6_1 || defined USE_COMPILER_PGI_6_2 || defined USE_COMPILER_PGI_7_0
    call getenv(trim(svar), svalueInEnv)
    if (trim(svalueInEnv) .eq. "") then
      nstatus = 1
    else
      nstatus = 0
    endif
#else
    call get_environment_variable(trim(svar), svalueInEnv, status=nstatus)
#endif

    select case (nstatus)
    case (0)
      ! Copy string only up to first whitespace character
!      read(svalueInEnv, '(A)') sresult

      ! Copy complete string
      sresult = svalueInEnv
      sys_getenv_string = .true.

    case (1)
      ! Environment variable does not exist
      sresult = ""
      sys_getenv_string = .false.

    case default
      !  2: Processor does not support environment variables
      ! >2: Some error occurred
      ! -1: variable svalueInEnv too short to absorb environment variables` content
      sresult = ""
      sys_getenv_string = .false.

    end select

  end function sys_getenv_string
  
! ****************************************************************************************

!<function>
  integer function sys_ncommandLineArgs()

  !<description>
    ! Calculates the number of command line arguments.
  !</description>
!</function>

#ifndef HAS_INTRINSIC_IARGC
    ! Definition of iargc needed for
    ! * Sun Fortran 95 8.1,
    ! * Compaq Fortran Compiler X5.4A-1684-46B5P,
    ! * Portland Group pgf90 6.0-5
    integer(I32) :: iargc
    external iargc ! generic Fortran routine to get the arguments from program call
#endif

    sys_ncommandLineArgs=iargc()

  end function
  
! ****************************************************************************************

!<subroutine>
  subroutine sys_getcommandLineArg(iarg,soption,svalue,iformat,sdefault)

  !<description>
    ! Fetches command line argument iarg from the command line.
    !
    ! The return value of this function depends on the format of the command line
    ! argument. There are three cases:
    !
    ! a) Simple option: "option", or svalue not specified
    !  Here, soption = "option" and
    !        svalue  = "".
    !        iformat = 0.
    !  Can be used to store e.g. paths like "./data".
    !
    ! b) Short options: "-option" or "-option=value".
    !  Here, soption = "option" and
    !        svalue  = "value" or "" if no value is specified.
    !        iformat = 1.
    !
    ! c) long options: "--option" or "--option=value".
    !  Here, soption = "option" and
    !        svalue  = "value" or "" if no value is specified.
    !        iformat = 2.
  !</description>
  
  !<input>
    ! Index of the command line argument. Must be in the range 1..sys_ncommandLineArgs().
    integer, intent(in) :: iarg

    ! OPTIONAL: A default value for the command line argument the iarg command
    ! line parameter does not exist.
    character(len=*), intent(in), optional :: sdefault
  !</input>
  
  !<output>
    ! The command line argument.
    character(len=*), intent(out) :: soption
    
    ! OPTIONAL: Value of the option
    character(len=*), intent(out), optional :: svalue

    ! OPTIONAL: Type of the command line argument.
    ! A value -1 indicates that the command line arguzment does not exist and no sdefault
    ! is specified.
    ! A value 0 indicates a direct option.
    ! A value 1 indicates that the command line parameter
    ! is of short form ("-key" or "-key=value").
    ! A value 2 indicates that the command line parameter
    ! is of long form ("--key" or "--key=value").
    integer, intent(out), optional :: iformat
  !</output>
  
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: stmp
    integer :: istmplen,idx

    if ((iarg .lt. 1) .or. (iarg .gt. sys_ncommandLineArgs())) then
      ! Return the default or an empty string.
      if (present(sdefault)) then
        stmp = sdefault
      else
        soption = ''
        if (present(iformat)) iformat = -1
        if (present(svalue)) svalue = ''
        return
      end if
    else
      ! Get the option -- at first to stmp.
      call getarg(iarg,stmp)
    end if

    istmplen = len(stmp)

    if (stmp(1:2).eq."--") then

      idx=3
      do while ((stmp(idx:idx).ne.'=').and.(idx.lt.istmplen))
        idx=idx+1
      enddo
  
      soption = stmp(3:idx-1)
      if (present(svalue)) svalue = stmp(idx+1:)
      if (present(iformat)) iformat = 2

    else if (stmp(1:1).eq."-") then

      idx = 2
      do while ((stmp(idx:idx).ne.'=').and.(idx .lt. istmplen))
        idx = idx+1
      enddo

      soption = stmp(2:idx-1)
      if (present(svalue)) svalue = stmp(idx+1:)
      if (present(iformat)) iformat = 1

    else

      soption = stmp
      if (present(svalue)) svalue = ""
      if (present(iformat)) iformat = 0

    endif

  end subroutine

! ****************************************************************************************

!<function>

  character (len=32) function sys_silsb(ivalue) result(soutput)

!<description>
    ! This routine converts an integer value to a string of length 32
    ! using LSB representation
!</description>

!<result>
    ! LSB Bit representation of the integer value with 32 bits.
!</result>

!<input>

    ! value to be converted
    integer, intent(in) :: ivalue

!</input>
!</function>

    ! local variables
    integer :: i
    
    do i = 1, min(32, bit_size(ivalue))
      if (btest(ivalue, i-1)) then
        soutput(i:i) = '1'
      else
        soutput(i:i) = '0'
      end if
    end do

    do i = min(32, bit_size(ivalue))+1, 32
      soutput(i:i) = '0'
    end do

  end function

! ****************************************************************************************

!<function>

  character (len=32) function sys_simsb(ivalue) result(soutput)

!<description>
    ! This routine converts an integer value to a string of length 32
    ! using MSB representation
!</description>

!<result>
    ! MSB Bit representation of the integer value with 32 bits.
!</result>

!<input>

    ! value to be converted
    integer, intent(in) :: ivalue

!</input>
!</function>

    ! local variables
    integer :: i
    
    do i = 1, min(32, bit_size(ivalue))
      if (btest(ivalue, i-1)) then
        soutput(32-i+1:32-i+1) = '1'
      else
        soutput(32-i+1:32-i+1) = '0'
      end if
    end do

    do i = min(32, bit_size(ivalue))+1, 32
      soutput(32-i+1:32-i+1) = '0'
    end do

  end function

! ****************************************************************************************

!<subroutine>

  subroutine sys_dequote (sstring)
  
!<description>
  ! Removes possible quotation marks around a string.
!</description>
  
!<inputoutput>
  ! String to de-quote
  character(len=*), intent(inout) :: sstring
!</inputoutput>

!</subroutine>
  
  character(len=len(sstring)+1) :: sstring2
  
    ! Adjust the string
    sstring2=trim(adjustl(sstring))
    
    ! Does the string start with a quotation mark?
    if ((sstring2(1:1) .eq. "'") .or. &
        (sstring2(1:1) .eq. """")) then
      ! Re-read the string, remove them.
      read (sstring2,*) sstring
    else
      ! Just transfer the string, it is ok.
      sstring = sstring2
    end if
  
  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine sys_stringToCharArray (sstring,schararray,slength)
  
!<description>
  ! Converts a string to a character array.
!</description>
  
!<input>
  ! String to convert
  character(len=*), intent(in) :: sstring
  
  ! OPTIONAL: Length of the string.
  ! If not present, the default string length is used.
  integer, intent(in), optional :: slength
!</input>

!<output>
  ! Character array that receives the converted string.
  ! Must be at least as long as the string or as slength.
  character, dimension(:), intent(out) :: schararray
!</output>

!</subroutine>
  
    integer :: i,j
    
    if (present(slength)) then
      j = slength
    else
      j = len_trim(sstring)
    end if
    
    ! Copy all characters.
    do i=1,j
      schararray(i) = sstring(i:i)
    end do
    
    ! Fill up the rest with spaces. This emulates a string copy.
    schararray(j+1:) = ' '
  
  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine sys_charArrayToString (schararray,sstring,slength)
  
!<description>
  ! Converts a character array to a string.
!</description>
  
!<input>
  ! Character array to convert
  character, dimension(:), intent(in) :: schararray

  ! OPTIONAL: Length of the string.
  ! If not present, the default string length is used.
  integer, intent(in), optional :: slength
!</input>

!<output>
  ! Character array that receives the converted string.
  ! Must be at least as long as the character array or slength.
  character(len=*), intent(out) :: sstring
!</output>

!</subroutine>
  
    integer :: i,j
    
    if (present(slength)) then
      j = slength
    else
      j = size(schararray)
    end if
    
    ! Copy all characters.
    sstring = ""
    do i=1,j
      sstring(i:i) = schararray(i)
    end do
  
  end subroutine

end module fsystem
