!##############################################################################
!# ****************************************************************************
!# <name> fsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains system routines like time measurement,          
!# string/value conversions, auxiliary routines and several sorting     
!# routines.
!#
!# On start of the main program, the routine system_init() must be called
!# once to initialise internal values!
!#
!# The following routines can be found here
!#
!#  1.) system_init
!#      -> Initialise system-wide settings
!#
!#  2.) sys_toupper
!#      -> Convert a string to uppercase
!#
!#  3.) sys_upcase
!#      -> Convert a string to uppercase, function version
!#
!#  4.) sys_throwFPE
!#      -> Throw a floating point exception
!#
!#
!#  5.) sys_getFreeUnit
!#      -> Determine a free file handle for use in an OPEN() command
!#
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
  
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(13,307)
  INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,63)

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

!<constantblock description="system signals">

  INTEGER, PARAMETER :: SIGILL  =  4
  INTEGER, PARAMETER :: SIGTRAP =  5
  INTEGER, PARAMETER :: SIGABRT =  6
  INTEGER, PARAMETER :: SIGEMT  =  7
  INTEGER, PARAMETER :: SIGFPE  =  8
  INTEGER, PARAMETER :: SIGBUS  = 10
  INTEGER, PARAMETER :: SIGSEGV = 11
  
!</constantblock>

!</constants>

CONTAINS

!************************************************************************

!<subroutine>

  SUBROUTINE sys_throwFPE()
  
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

  SUBROUTINE system_init()

    !<description>
    !This subroutine initialises some internal data structures.
    !</description>

!</subroutine>

    !set value of Pi = 3.14..

    SYS_PI=DASIN(1.0_DP)*2.0_DP

  END SUBROUTINE system_init

!************************************************************************

!<subroutine>
  
  SUBROUTINE sys_toupper (str) 

  !<description>
  
  ! Convert a string to upper case.
  
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
  
!******************************************************************************

!<function>
  FUNCTION sys_upcase(sinput) RESULT(soutput)
  
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
    integer              :: idigits
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
    integer              :: idigits
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
    integer             :: idigits
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
    integer             :: idigits
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
    integer                  :: idigits
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
    integer                  :: idigits
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
    integer              :: idigits
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
    integer              :: idigits
    !</input>
!</function>

    soutput = adjustl(sys_sdE(dvalue, idigits))
  end function sys_sdEL


!************************************************************************


!<function>
  character (len=32) function sys_siL(ivalue, idigits) result(soutput)

    !<description>
    ! This routine converts an integer value to a string of length idigits,
    ! filled up with white spaces.
    !</description>

    !<result>
    ! String representation of the value (left-aligned),
    ! fixed length of 32 characters
    !</result>

    !<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer  :: idigits
    !</input>
!</function>

    soutput = adjustl(sys_si(ivalue, idigits))
  end function sys_siL


!************************************************************************


!<function>
  character (len=32) function sys_si0L(ivalue, idigits) result(soutput)

    !<description>
    ! This routine converts an integer value to a string of length idigits,
    ! filled up with zeros.
    !</description>

    !<result>
    ! String representation of the value (left-aligned),
    ! fixed length of 32 characters
    !</result>

    !<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer             :: idigits
    !</input>
!</function>

    soutput = adjustl(sys_si0(ivalue, idigits))
  end function sys_si0L


!************************************************************************


!<function>
  character (len=32) function sys_sliL(ivalue, idigits) result(soutput)

    !<description>
    ! This routine converts a long integer value to a string of length idigits.
    !</description>

    !<result>
    ! String representation of the value (left-aligned),
    ! fixed length of 32 characters
    !</result>

    !<input>

    ! value to be converted
    integer(I64), intent(in) :: ivalue

    !number of decimals
    integer                  :: idigits
    !</input>
!</function>

    soutput = adjustl(sys_sli(ivalue, idigits))
  end function sys_sliL


!************************************************************************


!<function>
  character (len=32) function sys_sli0L(ivalue, idigits) result(soutput)

    !<description>
    ! This routine converts a long integer value to a string of length idigits.
    !</description>

    !<result>
    ! String representation of the value (left-aligned),
    ! fixed length of 32 characters
    !</result>

    !<input>

    ! value to be converted
    integer(I64), intent(in) :: ivalue

    !number of decimals
    integer                  :: idigits
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