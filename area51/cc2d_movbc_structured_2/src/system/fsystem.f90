!########################################################################
!# FINITE ELEMENT ANALYSIS & SOLUTION TOOLS  F E A S T  (Release 1.0)   #
!#                                                                      #
!# Authors: Ch.Becker,S.Kilian,S.Turek                                  #
!#          Institute of Applied Mathematics & Simulation               #
!#          University of Dortmund                                      #
!#          D-44227 DORTMUND                                            #
!#                                                                      #
!########################################################################
!#                                                                      #
!# <name> fsystem </name>                                               #
!#                                                                      #
!#                                                                      #
!# <purpose>                                                            #
!# This module contains system routines like time measurement,          #
!# string/value conversions, auxiliary routines and several sorting     #
!# routines.                                                            #
!# </purpose>                                                           #
!#                                                                      #
!########################################################################

!# Current version: $Id: fsystem.f90,v 1.113 2006/03/21 09:31:44 goeddeke4 Exp $

!<!--
#include <feastdefs.h>
! -->

module fsystem

  implicit none

!<constants>

!<constantblock description="constants for logical values">

  ! logical value 'true'
  integer, parameter :: YES = 0

  ! logical value 'false'
  integer, parameter :: NO = 1

!</constantblock>



#ifdef SINGLEPREC
!<constantblock description="kind values for floats">

  ! kind value for double precision
  integer, parameter :: DP = 4

  ! kind value for single precision
  integer, parameter :: SP = 4
!</constantblock>

#else
  integer, parameter :: DP = selected_real_kind(13,307)
  integer, parameter :: SP = selected_real_kind(6,63)
#endif

!<constantblock description="kind values for integers">

  ! kind value for 32Bit integer
  integer, parameter :: I32 = selected_int_kind(8)

  ! kind value for 64Bit integer
  integer, parameter :: I64 = selected_int_kind(10)

!</constantblock>


!<constantblock description="sort algorithms">

  ! heap sort (default); reliable allrounder
  integer, parameter :: SORT_HEAP   = 0

  ! quicksort (cutoff = 50), then insertsort
  integer, parameter :: SORT_QUICK  = 1

  ! insertsort; for small arrays
  integer, parameter :: SORT_INSERT = 2
!</constantblock>


!<constantblock description="system flags">

  !constant for a system beep
  character(len=1), parameter :: BEEP = achar(7)

  !constant for breaking line in a string (useful for
  !comm_sendmasterstring(...))
  character(len=1), parameter :: NEWLINE = achar(10)

  !standard length for strings in FEAST
  integer, parameter :: SYS_STRLEN = 256

  ! maximum number of additional variables in the master file
  integer, parameter :: SYS_MAXNADDVARS = 64

  ! maximum number of scarc solver definition files
  integer, parameter :: SYS_MAXNSCARC = 8

  ! maximum number of tokens per line in scarc solver definition file
  integer, parameter :: SYS_MAXNSCARCTOKENS = 20

  ! deprecated
  integer, parameter :: SYS_FEAT = 1

  !mathematical constant Pi
  real(DP)           :: SYS_PI

  !flag for debug mode
  integer            :: SYS_DEBUG

  !maximal values for real variables
  real(DP)           :: SYS_MAXREAL = huge(1.0_DP)

  !maximal values for integer variables
  integer            :: SYS_MAXINT = huge(1)

  ! increment value = 1
  integer, parameter :: INCX = 1

  ! control value for sys_deltatime
  integer, parameter :: SYS_TIMERSTART = 0

  ! control value for sys_deltatime
  integer, parameter :: SYS_TIMERSTOP = 1

  ! flag for appending data to a file (used in io)
  integer, parameter :: SYS_APPEND = 0

  ! flag for replacing a file  (used in io)
  integer, parameter :: SYS_REPLACE = 1

  ! flag for GMV ouput
  integer,parameter  :: SYS_GMV = 0

  ! flag for AVS output
  integer,parameter  :: SYS_AVS = 1

  ! file extension for GMV files
  character(len=6), parameter :: SYS_GMV_SUFFIX = ".gmv"

  ! file extension for AVS files
  character(len=6), parameter :: SYS_AVS_SUFFIX = ".inp"

  ! maximum number of global solution in one gmv file
  integer, parameter :: SYS_MAXGLOABLSOL = 6

!</constantblock>

!<constantblock description="system signals">
  integer, parameter :: SIGILL = 4
  integer, parameter :: SIGTRAP = 5
  integer, parameter :: SIGABRT = 6
  integer, parameter :: SIGEMT = 7
  integer, parameter :: SIGFPE = 8
  integer, parameter :: SIGBUS = 10
  integer, parameter :: SIGSEGV = 11
!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! global configuration type
  type t_sysconfig

    ! project id
    character(len=SYS_STRLEN) :: sprojectID

    ! project directory
    character(len=SYS_STRLEN) :: sprojectDir

    ! log directory
    character(len=SYS_STRLEN) :: slogDir

    ! output level
    integer :: coutputLevel

    ! name of feast grid file
    character(len=SYS_STRLEN) :: sgridFile


    character(len=SYS_STRLEN) :: scfg_gridfile_fbc
    character(len=SYS_STRLEN) :: scfg_gridfile_fpart
    character(len=SYS_STRLEN) :: scfg_gridfile_fmesh
    character(len=SYS_STRLEN) :: scfg_gridfile_fgeo
    character(len=SYS_STRLEN) :: scfg_gridfile_ffgeo
    character(len=SYS_STRLEN) :: sabsolutePath

    ! loadbalancing mode
    integer :: cloadBalanceMode

    ! file name for loadbalancing
    character(len=SYS_STRLEN) :: sloadBalanceFile

    ! anisotropic refinement mode
    integer :: canisoRefMode

    ! use constant matrix vector multiplication if possible
    integer :: cmatConstMul

    ! use fast matrix assembly if possible
    integer :: cfastMatAssembly

    ! calculate norm inside smoother
    integer :: ccalcSmoothNorm

    ! long output format
    integer :: clongOutputFormat

    ! number of scarc definitions
    integer :: nscarc

    ! scarc file names
    character(len=SYS_STRLEN), dimension(SYS_MAXNSCARC) :: sscarcFileNames

    ! data output file format as string
    character(len=SYS_STRLEN) :: soutputDataFormat

    ! data output file format as constant
    integer :: coutputDataFormat

    ! data output file name
    character(len=SYS_STRLEN) :: soutputDataFileName

    ! data output level: l .gt. 0 : perform advanced output on level l
    !                    l .lt. 0 : perform standard output (with duplicate vertices)
    !                               on level -l
    !                    l=0 : perform no output
    integer :: ioutputDataLevel

    ! flag if advanced output shall be used (advanced means: on macro boundaries there
    ! are no vertex duplicates anymore!)
    logical :: badvancedOutput

    ! use stored solution
    integer :: cuseStoredSolution

    ! use stored matrix
    integer :: cuseStoredMatrix

    ! default size of storage heaps
    integer :: idefaultStorageHeapSize

    ! default number of storage descriptors
    integer :: idefaultNStorageDescr

    ! maximum multigrid level
    integer :: imaxMGLevel

    ! sbblas mv implementation id
    integer :: isbblasMV

    ! sbblas prslowl implmentation id
    integer :: isbblasPRSLOWL

    ! sbblas mvtri implementation id
    integer :: isbblasMVTRI

    ! sbblas prsline implementation id
    integer :: isbblasPRSLINE

    ! index of additonal config variables
    integer :: iaddvaridx

    ! name of additional config parameters
    character(len=SYS_STRLEN), dimension(SYS_MAXNADDVARS) :: saddVarName

    ! value of additional config parameters
    character(len=SYS_STRLEN), dimension(SYS_MAXNADDVARS) :: saddVarValue

    ! autopartition
    integer :: cautoPartition

    ! autopartition with weighting
    integer :: cautoPartitionWeighted

    ! use extra storage for matrix edges
    integer :: cuseExtraStorageMatrixEdges
  end type t_sysconfig
!</typeblock>

!<typeblock>
  ! iteration descriptor
  type t_iterator
    ! start value
    integer :: istart

    ! stop value
    integer :: istop

    ! step value
    integer :: istep
  end type t_iterator
!</typeblock>


!<typeblock>
  ! simulation of an array of double pointers
  type t_realPointer
    real(DP), dimension(:), pointer :: ptr
  end type t_realPointer
!</typeblock>

!</types>


!<globals>

  real(DP) :: sys_dtimeMax

  ! global system configuration
  type (t_sysconfig) :: sys_sysconfig

  ! parallel mode 0 = single, 1=double binary (deprecated)
  integer :: sys_parmode

!</globals>


  real(DP) :: ars

!  external ztime

  contains

! The following function does not do what it should do. The string sresult
! does not really have the length len_trim(adjustl(sstring)) at the end, and
! so it contains trailing spaces again. There really seems to be NO WAY to
! return strings with variable length!

!   function sys_trimall(sstring) result(sresult)
!
!    !This function removes leading and trailing spaces.
!
!
!    !string
!    character(len=*) :: sstring
!
!    !string without leading and trailing spaces
!
!     character(len=len_trim(adjustl(sstring))) :: sresult
!
!     sresult=trim(adjustl(sstring))
!
!   end function


!************************************************************************


!<subroutine>
  subroutine sys_throwFPE()

    !<description>
    !This routine throws a floating point exception for debugging purposes
    !to prevent the debugger to exit the program.
    !</description>

!</subroutine>

    integer :: i1,i2

    i1=1
    i2=0

    i1=i1/i2

  end subroutine sys_throwFPE


!************************************************************************


!<function>
  character(len=SYS_STRLEN) function sys_sgetOutputFileBody()

    !<description>
    ! This routine provides the body of the outputfile name defined in the master.dat.
    ! Example : "solution.gmv", the routine provides "solution".
    !</description>

    !<result>
    !filename body
    !</result>

!</function>

    sys_sgetOutputFileBody = sys_sysconfig%soutputDataFileName

  end function sys_sgetOutputFileBody

!************************************************************************

!<function>
  character(len=6) function sys_sgetOutputFileSuffix()

    !<description>
    ! This routine provides the body of the outputfile name defined in the master.dat.
    ! Example : "solution.gmv", routine provides ".gmv".
    !</description>

    !<result>
    !filename extension
    !</result>

!</function>

    if (sys_sysconfig%coutputDataFormat .eq. SYS_GMV) then
      sys_sgetOutputFileSuffix = SYS_GMV_SUFFIX
    else if (sys_sysconfig%coutputDataFormat .eq. SYS_AVS) then
      sys_sgetOutputFileSuffix = SYS_AVS_SUFFIX
    endif

  end function sys_sgetOutputFileSuffix

!************************************************************************

!<function>
  character(len=SYS_STRLEN) function sys_sgetOutputFileName()

    !<description>
    ! This routine composes the complete name of the visualisation output file
    ! from basename specified in a configuration file (master.dat) and the
    ! appropriate suffix based on the requested visualisation output format.
    !</description>

    !<result>
    !filename
    !</result>

!</function>

    sys_sgetOutputFileName = trim(sys_sysconfig%soutputDataFileName) // &
                             trim(sys_sgetOutputFileSuffix())

  end function sys_sgetOutputFileName


!************************************************************************************


!<subroutine>
  subroutine sys_version(ifeastVersionHigh, ifeastVersionMiddle, ifeastVersionLow, &
                         sreldate)

    !<description>
    !This subroutine returns the library version information.
    !</description>

    !<output>

    ! high version number
    integer :: ifeastVersionHigh

    ! middle version number
    integer :: ifeastVersionMiddle

    ! low version number
    integer :: ifeastVersionLow

    ! release date
    character(len=*) :: sreldate

    !</output>

!</subroutine>

    ifeastVersionHigh=0
    ifeastVersionMiddle=9
    ifeastVersionLow=0

    sreldate="22.11.2005 RC1"

  end subroutine sys_version


!************************************************************************************


!<subroutine>
  subroutine sys_getBuildCPU(scpu)

    !<description>
    !This subroutine returns informations about the built cpu.
    !</description>

    !<output>

    ! build architecture
    character(len=SYS_STRLEN) :: sarch

    ! build cpu
    character(len=*) :: scpu

    ! build operating system
    character(len=SYS_STRLEN) :: sos

    ! use MPI environment
    character(len=SYS_STRLEN) :: smpienv

    ! build compiler
    character(len=SYS_STRLEN) :: scompiler

    ! used BLAS
    character(len=SYS_STRLEN) :: sblas

    !</output>

!</subroutine>

    include "buildconf.h"

  end subroutine sys_getBuildCPU


!************************************************************************************


!<subroutine>
  subroutine sys_getBuildCompiler(scompiler)

    !<description>
    !This subroutine returns informations about the build architecture.
    !</description>

    !<output>

    ! build architecture
    character(len=SYS_STRLEN) :: sarch

    ! build cpu
    character(len=SYS_STRLEN) :: scpu

    ! build operating system
    character(len=SYS_STRLEN) :: sos

    ! use MPI environment
    character(len=SYS_STRLEN) :: smpienv

    ! build compiler
    character(len=*) :: scompiler

    ! used BLAS
    character(len=SYS_STRLEN) :: sblas

    !</output>

!</subroutine>

    include "buildconf.h"

  end subroutine sys_getBuildCompiler


!************************************************************************************


!<subroutine>
  subroutine sys_getBuildArch(sarch)

    !<description>
    !This subroutine returns informations about the build architecture.
    !</description>

    !<output>

    ! build architecture
    character(len=*) :: sarch

    ! build cpu
    character(len=SYS_STRLEN) :: scpu

    ! build operating system
    character(len=SYS_STRLEN) :: sos

    ! use MPI environment
    character(len=SYS_STRLEN) :: smpienv

    ! build compiler
    character(len=SYS_STRLEN) :: scompiler

    ! used BLAS
    character(len=SYS_STRLEN) :: sblas

    !</output>

!</subroutine>

    include "buildconf.h"

  end subroutine sys_getBuildArch


!************************************************************************************


!<subroutine>
  subroutine sys_getBuildEnv(sarch, scpu, sos, smpienv, scompiler, sblas)

    !<description>
    !This subroutine returns informations about the built environment.
    !</description>

    !<output>

    ! build architecture
    character(len=*) :: sarch

    ! build cpu
    character(len=*) :: scpu

    ! build operating system
    character(len=*) :: sos

    ! use MPI environment
    character(len=*) :: smpienv

    ! build compiler
    character(len=*) :: scompiler

    ! used BLAS
    character(len=*) :: sblas

    !</output>

!</subroutine>

    include "buildconf.h"

  end subroutine sys_getBuildEnv


!************************************************************************


!<subroutine>
  subroutine sys_init()

    !<description>
    !This subroutine initialises some internal data structures.
    !</description>

    !<global variable="sys_dtimeMax, SYS_MATCONST, SYS_FASTASS, SYS_CALCSMOOTHNORM, SYS_PI">
    !</global>
!</subroutine>

    integer :: icount ! current system time
    integer :: irate  ! approx. number of system clock ticks per second
    integer :: icmax  ! largest possible value of icount

    !system_clock is not a FEAST, but a basic FORTRAN 90 routine
    call system_clock(icount,irate,icmax)

    !maximal measurable time span in seconds (system-dependend)
    sys_dtimeMax = real(icmax,DP)/real(irate,DP)

    !these values will be overwritten in master_start in mastermod.f90
    !by the values defined in master.dat
    sys_sysconfig%cmatConstMul= NO
    sys_sysconfig%cfastMatAssembly = NO
    sys_sysconfig%ccalcSmoothNorm = NO

    !set value of Pi = 3.14..
#ifdef SINGLEPREC
    SYS_PI=asin(1.0_DP)*2.0_DP
#else
    SYS_PI=dasin(1.0_DP)*2.0_DP
#endif

    if (sys_getenv_int("FEASTDEBUGSESSION",icount)) then
      SYS_DEBUG = YES
    else
      SYS_DEBUG = NO
    endif

    sys_sysconfig%cuseExtraStorageMatrixEdges=NO

  end subroutine sys_init


!************************************************************************


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


!************************************************************************


!<subroutine>
  function sys_inumberOfDigits(inumber)
    !<description>
    !This function determines the number of digits a given number has.
    !</description>

    !<input>
    ! name of the master file
    integer :: inumber
    !</input>

    !<output>

    !the number of digits of the given number
    integer :: sys_inumberOfDigits
    !</output>
!</subroutine>

    sys_inumberOfDigits = log(real(inumber, SP)) / log(10.0_SP) + 1

  end function sys_inumberOfDigits


!****************************************************************************


!<subroutine>
  subroutine sys_getMasterFile(smaster)

    !<description>
    !This subroutine returns the name of the master .dat file.
    !
    !The default dat-file is "master.dat". If a different file shall
    !be used as dat-file, its name has to be appended to the program call:\\
    !
    !Example: mpirun -np x program-name [dat-file]\\
    !
    !              x - number of processors to be used\\
    !
    !       dat-file - name of the dat-file to be used instead of master.dat
    !</description>

    !<output>

    ! name of the master file
    character (len=*) :: smaster

    !</output>

!</subroutine>

#ifndef HAS_INTRINSIC_IARGC
    !Definition of iargc needed for
    !* Sun Fortran 95 8.1,
    !* Compaq Fortran Compiler X5.4A-1684-46B5P,
    !* Portland Group pgf90 6.0-5
    integer :: iargc
    external iargc !generic Fortran routine to get the arguments from program call
#endif

    if (iargc() .ne. 1) then  !if no optional argument
      smaster="master.dat"  !default name
    else
      call getarg(1,smaster)!else get name of master file
    endif

  end subroutine sys_getMasterFile


!****************************************************************************


!<subroutine>
  recursive subroutine sys_tokeniser(sbuffer, nmaxtokcount, stokens, itokcount)

    !<description>
    ! This routine splits a string into substrings, divided by blank,
    ! comma or equation sign.
    !</description>

    !<input>

    ! string to be scanned
    character(len=*) :: sbuffer

    ! maximal number of substrings
    integer :: nmaxtokcount

    !</input>

    !<output>
    ! array of substrings
    character(len=*), dimension(:) :: stokens

    ! number of tokens
    integer, intent(out) :: itokCount
    !</output>
!</subroutine>

    !Copy of input string plus additional blank to have a
    !token delimiter even after last token.
    character(len=len(sbuffer) + 1) :: sbufCopy

    !loop counter
    integer :: i, idx

    ! Position of first non-blank in string buffer
    integer :: itokensStart

    ! Position of last non-blank in string buffer
    integer :: itokensEnd

    ! Start position of current token
    integer :: icurrTokStart

    !Create trimmed copy of string buffer
    sbufCopy = trim(adjustl(sbuffer))

    itokCount = 0
    itokensStart  = 1
    itokensEnd    = len(trim(sbufCopy)) + 1  ! +1 to find also a delimiter
                                             ! after the last token
    icurrTokStart = itokensStart
    do idx = itokensStart, itokensEnd
      ! Token delimiter reached? Then save it.
      if ( (sbufCopy(idx:idx) .eq. '=') .or. &
           (sbufCopy(idx:idx) .eq. '@') .or. &
           (sbufCopy(idx:idx) .eq. ',') .or. &
           (sbufCopy(idx:idx) .eq. ' ') ) then
        stokens(itokCount+1) = sbufCopy(icurrTokStart:idx - 1)
        icurrTokStart = idx + 1
        itokCount = itokCount + 1
      endif

      if (itokCount .eq. nmaxtokcount) then
        return
      endif
    enddo

  end subroutine sys_tokeniser


!************************************************************************


!<function>
  logical function sys_getenv_int(svar,ivalue)

    !<description>
    ! This functions returns the integer value of a given
    ! enviroment variable. The routine returns .TRUE., if the variable
    ! exists, otherwise .FALSE. .
    !</description>

    !<result>
    ! exit status
    !</result>

    !<input>

    ! name of the enviroment variable
    character(len=*) :: svar

    !</input>

    !<output>

    ! value of the enviroment variable
    integer :: ivalue
    !</output>
!</function>

    character(len=128) :: svalue !output string

    call getenv_wrap(trim(svar),svalue,len(svar),128)

    if (trim(svalue).eq."") then
      ivalue=0
      sys_getenv_int=.FALSE.
    else
      read(svalue,*) ivalue
      sys_getenv_int=.TRUE.
    endif

  end function sys_getenv_int


!************************************************************************


!<function>
  logical function sys_getenv_real(svar,dvalue)

    !<description>
    ! This functions returns the real value of a given
    ! enviroment variable. The routine returns .TRUE., if the variable
    ! exists, otherwise .FALSE. .
    !</description>

    !<result>
    ! exit status
    !</result>

    !<input>

    ! name of the enviroment variable
    character(len=*) :: svar

    !</input>

    !<output>
    real(DP) :: dvalue
    !</output>
!</function>

    character(len=128) :: svalue

    call getenv_wrap(trim(svar),svalue,len(svar),128)

    if (trim(svalue).eq."") then
      dvalue=0.0
      sys_getenv_real=.FALSE.
    else
      read(svalue,*) dvalue
      sys_getenv_real=.TRUE.
    endif

  end function sys_getenv_real


!************************************************************************


!<function>
  logical function sys_getenv_string(svar,sresult)

    !<description>
    ! This functions returns the string value of a given
    ! enviroment variable. The routine returns .TRUE., if the variable
    ! exists, otherwise .FALSE. .
    !</description>

    !<result>
    ! exit status
    !</result>

    !<input>

    ! name of the enviroment variable
    character(len=*) :: svar

    !</input>

    !<output>

    character(len=*) :: sresult

    !</output>
!</function>

    character(len=SYS_STRLEN) :: svalueInEnv

    call getenv_wrap(trim(svar), svalueInEnv, len(svar), SYS_STRLEN)

    if (trim(svalueInEnv).eq."") then
      sresult=""
      sys_getenv_string = .FALSE.
    else
      read(svalueInEnv,'(A)') sresult
      sys_getenv_string = .TRUE.
    endif

  end function sys_getenv_string


!************************************************************************


!<function>
  integer function sys_readpar_int(svalue)

    !<description>
    ! This routine checks, if the given string starts with an \$-character
    ! and returns in this case the integer value of the enviroment variable.
    ! Otherwise it returns the integer value of the string itself.
    !</description>

    !<result>
    ! integer value of the string or environment variable
    !</result>

    !<input>

    ! string with value or variable name
    character(len=*) :: svalue

    !</input>
!</function>

    integer :: itmp

    logical :: bexists

    if (svalue(1:1).eq.'$') then
      bexists = sys_getenv_int(svalue(2:),itmp)

      if (.not.bexists) then
        write (6,*) "*** WARNING! env variable ",trim(svalue(2:))," not set. ***"
      endif

    else
      itmp = sys_stringToInt(svalue)
    endif

    sys_readpar_int=itmp

  end function sys_readpar_int


!************************************************************************


!<function>
 real(DP) function sys_readpar_real(svalue)

    !<description>
    ! This routine checks, if the given string starts with an \$-character
    ! and returns in this case the integer value of the enviroment variable.
    ! Otherwise it returns the real value of the string itself.
    !</description>

    !<result>
    ! integer value of the string or environment variable
    !</result>

    !<input>

    ! string with value or variable name
    character(len=*) :: svalue
    !</input>

!</function>

    real(DP) :: dtmp

    logical :: bexists

    if (svalue(1:1).eq.'$') then
      bexists = sys_getenv_real(svalue(2:),dtmp)
      if (.not.bexists) then
        write (6,*) "*** WARNING! env variable ",trim(svalue(2:))," not set. ***"
      endif
    else
      dtmp = sys_stringToReal(svalue)
    endif

    sys_readpar_real=dtmp

  end function sys_readpar_real


!************************************************************************


!<function>
  character(len=256) function sys_readpar_string(svalue)

    !<description>
    ! This routine checks, if the given string starts with an \$-character
    ! and returns in this case the string value of the enviroment variable.
    ! Otherwise it returns the string value of the string itself.
    !</description>

    !<result>
    ! integer value of the string or environment variable
    !</result>

    !<input>

    ! string with value or variable name
    character(len=*) :: svalue

    !</input>
!</function>

    character(len=256) :: stmp

    logical :: bexists

    if (svalue(1:1).eq.'$') then
      bexists = sys_getenv_string(svalue(2:),stmp)
      if (.not.bexists) then
        write (6,*) "*** WARNING! env variable ",trim(svalue(2:))," not set. ***"
      endif
    else
      read(svalue,'(A)') stmp
    endif

    sys_readpar_string=stmp

  end function sys_readpar_string


!************************************************************************


!<function>
  function sys_upcase(sinput) result(soutput)
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
    integer :: i

    soutput = " "   !initialise string
    do i = 1,len(sinput)
       if(sinput(i:i) .ge. "a" .and. sinput(i:i) .le. "z") then
          soutput(i:i) = achar(iachar(sinput(i:i)) - 32)
       else
          soutput(i:i) = sinput(i:i)
       end if
    end do
  end function sys_upcase


!************************************************************************


!<function>
  integer function sys_stringToInt(svalue)

    !<description>
    ! This routine tries to convert a string to an integer value. If the conversion fails,
    ! the return value is set to SYS_MAXINT.
    !</description>

    !<input>

    !string to be converted
    character(*), intent(in) :: svalue

    !</input>

    !<result>
    !resulting value ( = SYS_MAXINT if conversion fails)
    !</result>
!</function>

    integer :: istatus !status variable indicating if error occurs

    !read svalue and try to write into iValue
    read(svalue, fmt=*, iostat=istatus) sys_stringToInt

    if (istatus .ne. 0) then
      sys_stringToInt = SYS_MAXINT
      write (6,*) "*** WARNING! String '" // trim(svalue) // &
                  "' could not be converted to integer. ***"
    endif

  end function sys_stringToInt


!************************************************************************


!<function>
  real(DP) function sys_stringToReal(svalue)

    !<description>
    ! This routine tries to convert a string to a real value. If the conversion fails,
    ! the return value is set to SYS_MAXREAL.
    !</description>

    !<input>

    !string to be converted
    character(*), intent(in) :: svalue

    !</input>

    !<result>
    !resulting value ( = SYS_MAXREAL if conversion fails)
    !</result>
!</function>

    integer :: istatus !status variable indicating if error occurs

    !read svalue and try to write into dValue
    read(svalue, fmt=*, iostat=istatus) sys_stringToReal

    if (istatus .ne. 0) then
      sys_stringToReal = SYS_MAXREAL
      write (6,*) "*** WARNING! String '" // trim(svalue) // &
                  "' could not be converted to real. ***"
    endif

  end function sys_stringToReal


!************************************************************************


!******************************************************************************
! Routines for file handling
!******************************************************************************


!<subroutine>
  subroutine sys_getNextEntry(iunit, sbuffer)

    !<description>
    !This routine reads an item from the file connected to unit iunit and writes it
    !into the buffer sbuffer. Items are assumed to be separated by spaces, comma, tabs.
    !Anything from hash character \# aka Lattenkreuz till the EOL is ignored (assuming
    !fortrans eor=EOL). The escape symbol is the backslash.\\
    !Author: Jaroslav
    !</description>

    !<input>

    !unit connected to the file to read from
    integer, intent(in) :: iunit
    !</input>

    !<output>

    !output string
    character(len=*), intent(out) :: sbuffer
    !</output>
!</subroutine>

    character(len=1) :: saux

    integer :: i !index variable

    ! variables to expand environment variables on-the-fly when found
    character(len=SYS_STRLEN), save :: sauxEnv
    logical :: bfoundInEnv

    sbuffer = ""
    i=1

    do
      !read character
      !If EOF is reached then go to line 99. If end of line (eor = end of record)
      !is reached goto 10, which means that empty lines are skipped.
10    read(unit=iunit, fmt='(a)', advance='no', end=99, eor=10) saux

      !exit loop if charac is not a separating character
      if((saux .ne. ' ') .and. (saux .ne. ',') .and. (saux .ne. achar(9)) .and. &
         (saux .ne. '#')) then
        exit
      endif

      !If # is found, read until end of line.
      if(saux .eq. '#') then
        do
          read(unit=iunit, fmt='(a)', advance='no', end=99, eor=10) saux
        enddo
      endif
    enddo

    !Read the remaining characters of the item
    do
      sbuffer(i:i)=saux; i=i+1
      read(unit=iunit, fmt='(a)', advance='no', end=99, eor=99) saux

      ! achar(9)=tabulator
      if((saux .eq. ' ') .or. (saux .eq. ',') .or. (saux .eq. achar(9))) then
        exit
      endif
      if(saux.eq.'#') then
        write(*,*) sbuffer(i-1:i-1)

        !'#' marks comments
        !unless the previous character is the escape character (achar(92)=backslash)
        if ((i .gt. 1) .and. (sbuffer(i-1:i-1) .ne. achar(92))) then
          do
            read(unit=iunit, fmt='(a)', advance='no', end=99, eor=99) saux
          enddo
          write(*,*) "Error in routine sys_getNextEntry()!"
        else

          !to overwrite escape character
          i=i-1
        endif
      endif
    enddo
99  sbuffer(i:i)=' '

    !Check whether item is an environment variable
    if (sbuffer(1:1) .eq. '$') then
      bfoundInEnv = .FALSE.
      ! Retrieve value of environment variable
      ! (Do not forget to cut the dollar sign.)
      bfoundInEnv = sys_getenv_string(trim(sbuffer(2:)), sauxEnv)
      if (bfoundInEnv) then
        ! Do not return the environment variable, return its expanded value!
        sbuffer = trim(sauxEnv)
      else
        write (6,*) "*** WARNING! Could not expand environment variable '" // &
                    trim(sbuffer) // "'! ***"
      endif
    endif

  end subroutine sys_getNextEntry


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
    do itry = 1,10000
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
    integer :: iunit

    !</input>
!</function>

#ifdef HAS_FLUSH
    call flush(iunit)
#endif

  end subroutine sys_flush


!************************************************************************


!******************************************************************************
! Sorting routines
!******************************************************************************


!<subroutine>
  subroutine sys_int_sort(Iarray, csortMethod, Imapping)

    !<description>
    !Sorting routine for integer data type. If more than one vector must be
    !sorted, or if several vectors may be sorted the same way, the mapping
    !array Imapping can be computed. This array must be created outside the routine
    !and must have the same length as Iarray. Then, the other vectors are computed
    !by Vec(Imapping(i)). The array Imapping is an optional argument. The sorting
    !method used is told by the parameter csortMethod. If this optional argument
    !is not used, the routine performs heapsort.
    !</description>

    !<input>

    ! sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer,optional::csortMethod

    !</input>

    !<inoutput>

    ! integer array to be sorted
    integer,dimension(:) :: Iarray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer,dimension(:), optional:: Imapping
    !</inoutput>
!</subroutine>


    !if the optional argument csortMethod is present
    if (present(Imapping)) then
      if(present(csortMethod)) then
        select case (csortMethod)
        case(SORT_HEAP)
          call heapsort(Iarray, Imapping)
        case(SORT_QUICK)
          call random_seed
          call quicksort(Iarray, Imapping)
          call insertsort(Iarray, Imapping)
        case(SORT_INSERT)
          call insertsort(Iarray, Imapping)
        case default
          stop 'sys_i32_sort: unknown method'
        end select
      else
        call heapsort(Iarray, Imapping)
      endif
    else
      if(present(csortMethod)) then
        select case (csortMethod)
        case(SORT_HEAP)
          call heapsort(Iarray)
        case(SORT_QUICK)
          call random_seed
          call quicksort(Iarray)
          call insertsort(Iarray)
        case(SORT_INSERT)
          call insertsort(Iarray)
        case default
          stop 'sys_i32_sort: unknown method'
        end select
      else
        call heapsort(Iarray)
      endif
    endif
  contains

    !************************************************************************

    subroutine reheap(Iarray, istart, istop, Imapping)
      integer::istart,istop
      integer,dimension(1:istop)::Iarray
      integer,dimension(1:istop), optional:: Imapping
      integer::t1,t2
      integer :: i1, i2
      integer::i,j
      !if(istop.eq.start) return ! nothing to correct

      if (present(Imapping)) then

        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j>0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        end do
        
        !search correct position along the path
        do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
          i=ishft(i,-1)
        end do
        
        !move data
        t1=Iarray(istart)
        i1=Imapping(istart)
        
        do while (i>=istart)
          t2=Iarray(i)
          Iarray(i)=t1
          t1=t2
          
          i2=Imapping(i)
          Imapping(i)=i1

          i1=i2
          i=ishft(i,-1)
        end do
      else
        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j>0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        end do
        
        !search correct position along the path
        do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
          i=ishft(i,-1)
        end do
        
        !move data
        t1=Iarray(istart)
        
        do while (i>=istart)
          t2=Iarray(i)
          Iarray(i)=t1
          t1=t2
          i=ishft(i,-1)
        end do
      endif

    end subroutine reheap

    !************************************************************************

    subroutine heapsort(Iarray, Imapping)
      integer,dimension(:)::Iarray
      integer,dimension(:), optional:: Imapping
      integer::t
      integer :: t2
      integer :: i,n
      n = ubound(Iarray,1)
      ! heap creation phase (Maxheap)

      if (present(Imapping)) then
        do i=ishft(n,-1),1,-1
          call reheap(Iarray,i,n,Imapping)
        end do
        ! selection phase
        do i=n,2,-1
          t2=Imapping(i)
          Imapping(i)=Imapping(1)
          Imapping(1)=t2

          t=Iarray(i)
          Iarray(i)=Iarray(1)
          Iarray(1)=t
          call reheap(Iarray,1,i-1,Imapping)
        end do
      else
        do i=ishft(n,-1),1,-1
          call reheap(Iarray,i,n)
        end do
        ! selection phase
        do i=n,2,-1
          
          t=Iarray(i)
          Iarray(i)=Iarray(1)
          Iarray(1)=t
          call reheap(Iarray,1,i-1)
        end do
      endif
    end subroutine heapsort

    !************************************************************************

    recursive subroutine quicksort(Iarray, Imapping)
      integer,dimension(:)::Iarray
      integer,dimension(:), optional:: Imapping
      integer::t,temp

      integer :: t2, temp2

      integer::l,u,i,j
      real:: r
      integer,parameter::cutoff=50
      l=1
      u=ubound(Iarray,1)

      if (present(Imapping)) then
        do while ((u-l)>=cutoff)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Iarray(i)
          Iarray(i)=Iarray(l)
          Iarray(l)=t
          
          t2=Imapping(i)
          Imapping(i)=Imapping(l)
          Imapping(l)=t2

          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Iarray(i).lt.t))
              i=i+1
            end do
            do while (Iarray(j).gt.t)
              j=j-1
            end do
            if (i.gt.j) exit
            temp=Iarray(i)
            temp2=Imapping(i)

            Iarray(i)=Iarray(j)
            Iarray(j)=temp

            Imapping(i)=Imapping(j)
            Imapping(j)=temp2

            j=j-1
          end do

          Iarray(l)=Iarray(j)
          Iarray(j)=t
          
          Imapping(l)=Imapping(j)
          Imapping(j)=t2
          
          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Iarray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            call quicksort(Iarray(l:j-1), Imapping(l:j-1))
            l=j+1
          endif
        end do
      else
        
        do while ((u-l)>=cutoff)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Iarray(i)
          Iarray(i)=Iarray(l)
          Iarray(l)=t
          
          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Iarray(i).lt.t))
              i=i+1
            end do
            do while (Iarray(j).gt.t)
              j=j-1
            end do
            if (i.gt.j) exit
            temp=Iarray(i)
            Iarray(i)=Iarray(j)
            Iarray(j)=temp
            j=j-1
          end do
          Iarray(l)=Iarray(j)
          Iarray(j)=t
          
          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Iarray(j+1:u))
            u=j-1
          else
            call quicksort(Iarray(l:j-1))
            l=j+1
          endif
        end do
        
      endif
    end subroutine quicksort

    !************************************************************************

    subroutine insertsort(Iarray, Imapping)
      integer,dimension(:)::Iarray
      integer,dimension(:), optional:: Imapping
      integer::t
      integer :: t2
      integer::i,j

      if (present(Imapping)) then
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          t2 = Imapping(i)
          j=i-1

          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          end do
          Iarray(j+2:i)=Iarray(j+1:i-1)
          Iarray(j+1)=t

          Imapping(j+2:i)=Imapping(j+1:i-1)
          Imapping(j+1)=t2
        end do
      else
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          j=i-1
          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          end do
          Iarray(j+2:i)=Iarray(j+1:i-1)
          Iarray(j+1)=t
        end do
      endif
    end subroutine insertsort

  end subroutine sys_int_sort


!************************************************************************


!<subroutine>
  subroutine sys_i32_sort(iarray, csortMethod, Imapping)


    !<description>
    !Sorting routine for i32 data type. If more than one vector must be
    !sorted, or if several vectors may be sorted the same way, the mapping
    !array Imapping can be computed. This array must be created outside the routine
    !and must have the same length as Iarray. Then, the other vectors are computed
    !by Vec(Imapping(i)). The array Imapping is an optional argument. The sorting
    !method used is told by the parameter csortMethod. If this optional argument
    !is not used, the routine performs heapsort.
    !</description>


    !<input>

    ! sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer,optional::csortMethod

    !</input>

    !<inoutput>

    ! integer array to be sorted
    integer(i32),dimension(:) :: Iarray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer,dimension(:), optional:: Imapping
    !</inoutput>
!</subroutine>

    !if the optional argument csortMethod is present
    if (present(Imapping)) then
      if(present(csortMethod)) then
        select case (csortMethod)
        case(SORT_HEAP)
          call heapsort(Iarray, Imapping)
        case(SORT_QUICK)
          call random_seed
          call quicksort(Iarray, Imapping)
          call insertsort(Iarray, Imapping)
        case(SORT_INSERT)
          call insertsort(Iarray, Imapping)
        case default
          stop 'sys_i32_sort: unknown method'
        end select
      else
        call heapsort(Iarray, Imapping)
      endif
    else
      if(present(csortMethod)) then
        select case (csortMethod)
        case(SORT_HEAP)
          call heapsort(Iarray)
        case(SORT_QUICK)
          call random_seed
          call quicksort(Iarray)
          call insertsort(Iarray)
        case(SORT_INSERT)
          call insertsort(Iarray)
        case default
          stop 'sys_i32_sort: unknown method'
        end select
      else
        call heapsort(Iarray)
      endif
    endif
  contains

    !************************************************************************

    subroutine reheap(Iarray, istart, istop, Imapping)
      integer::istart,istop
      integer(i32),dimension(1:istop)::Iarray
      integer,dimension(1:istop), optional:: Imapping
      integer(i32)::t1,t2
      integer :: i1, i2
      integer::i,j
      !if(istop.eq.start) return ! nothing to correct

      if (present(Imapping)) then

        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j>0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        end do
        
        !search correct position along the path
        do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
          i=ishft(i,-1)
        end do
        
        !move data
        t1=Iarray(istart)
        i1=Imapping(istart)
        
        do while (i>=istart)
          t2=Iarray(i)
          Iarray(i)=t1
          t1=t2
          
          i2=Imapping(i)
          Imapping(i)=i1
          i1=i2
          i=ishft(i,-1)
        end do
      else
        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j>0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        end do
        
        !search correct position along the path
        do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
          i=ishft(i,-1)
        end do
        
        !move data
        t1=Iarray(istart)
        
        do while (i>=istart)
          t2=Iarray(i)
          Iarray(i)=t1
          t1=t2
          i=ishft(i,-1)
        end do
      endif

    end subroutine reheap

    !************************************************************************

    subroutine heapsort(Iarray, Imapping)
      integer(i32),dimension(:)::Iarray
      integer,dimension(:), optional:: Imapping
      integer(i32)::t
      integer :: t2
      integer :: i,n
      n=ubound(Iarray,1)
      ! heap creation phase (Maxheap)

      if (present(Imapping)) then
        do i=ishft(n,-1),1,-1
          call reheap(Iarray,i,n,Imapping)
        end do
        ! selection phase
        do i=n,2,-1
          t2=Imapping(i)
          Imapping(i)=Imapping(1)
          Imapping(1)=t2

          t=Iarray(i)
          Iarray(i)=Iarray(1)
          Iarray(1)=t
          call reheap(Iarray,1,i-1,Imapping)
        end do
      else
        do i=ishft(n,-1),1,-1
          call reheap(Iarray,i,n)
        end do
        ! selection phase
        do i=n,2,-1
          
          t=Iarray(i)
          Iarray(i)=Iarray(1)
          Iarray(1)=t
          call reheap(Iarray,1,i-1)
        end do
      endif
    end subroutine heapsort

    !************************************************************************

    recursive subroutine quicksort(Iarray, Imapping)
      integer(i32),dimension(:)::Iarray
      integer,dimension(:), optional:: Imapping
      integer(i32)::t,temp

      integer :: t2, temp2

      integer::l,u,i,j
      real:: r
      integer,parameter::cutoff=50
      l=1
      u=ubound(Iarray,1)

      if (present(Imapping)) then
        do while ((u-l)>=cutoff)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Iarray(i)
          Iarray(i)=Iarray(l)
          Iarray(l)=t
          
          t2=Imapping(i)
          Imapping(i)=Imapping(l)
          Imapping(l)=t2

          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Iarray(i).lt.t))
              i=i+1
            end do
            do while (Iarray(j).gt.t)
              j=j-1
            end do
            if (i.gt.j) exit
            temp=Iarray(i)
            temp2=Imapping(i)

            Iarray(i)=Iarray(j)
            Iarray(j)=temp

            Imapping(i)=Imapping(j)
            Imapping(j)=temp2

            j=j-1
          end do

          Iarray(l)=Iarray(j)
          Iarray(j)=t
          
          Imapping(l)=Imapping(j)
          Imapping(j)=t2
          
          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Iarray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            call quicksort(Iarray(l:j-1), Imapping(j+1:u))
            l=j+1
          endif
        end do
      else
        
        do while ((u-l)>=cutoff)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Iarray(i)
          Iarray(i)=Iarray(l)
          Iarray(l)=t
          
          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Iarray(i).lt.t))
              i=i+1
            end do
            do while (Iarray(j).gt.t)
              j=j-1
            end do
            if (i.gt.j) exit
            temp=Iarray(i)
            Iarray(i)=Iarray(j)
            Iarray(j)=temp
            j= j-1
          end do
          Iarray(l)=Iarray(j)
          Iarray(j)=t
          
          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Iarray(j+1:u))
            u=j-1
          else
            call quicksort(Iarray(l:j-1))
            l=j+1
          endif
        end do
        
      endif
    end subroutine quicksort

    !************************************************************************

    subroutine insertsort(Iarray, Imapping)
      integer(i32),dimension(:)::Iarray
      integer,dimension(:), optional:: Imapping
      integer(i32)::t
      integer :: t2
      integer::i,j

      if (present(Imapping)) then
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          t2 = Imapping(i)
          j=i-1
          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          end do
          Iarray(j+2:i)=Iarray(j+1:i-1)
          Iarray(j+1)=t

          Imapping(j+2:i)=Imapping(j+1:i-1)
          Imapping(j+1)=t2
        end do
      else
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          j=i-1
          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          end do
          Iarray(j+2:i)=Iarray(j+1:i-1)
          Iarray(j+1)=t
        end do
      endif
    end subroutine insertsort

  end subroutine sys_i32_sort


!************************************************************************


!<subroutine>
  subroutine sys_i64_sort(Iarray, csortMethod)

    !<description>
    ! sort routine for integer64 type
    !</description>

    !<input>

    ! sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer,optional::csortMethod

    !</input>

    !<inoutput>

    ! integer array to be sorted
    integer(i64),dimension(:)::Iarray
    !</inoutput>
!</subroutine>

    !if the optional argument csortMethod is present
    if(present(csortMethod)) then
       select case (csortMethod)
       case(SORT_HEAP)
          call heapsort(Iarray)
       case(SORT_QUICK)
          call random_seed
          call quicksort(Iarray)
          call insertsort(Iarray)
       case(SORT_INSERT)
          call insertsort(Iarray)
       case default
          stop 'sys_i64_sort: unknown method'
       end select
    else
       call heapsort(Iarray)
    endif
  contains

    !************************************************************************

    subroutine reheap(Iarray, istart, istop)
      integer::istart,istop
      integer(i64),dimension(1:istop)::Iarray
      integer(i64)::t1,t2
      integer::i,j
      !if(istop.eq.istart) return ! nothing to correct

      !trace the path of the bigger children
      i=istart
      j=ishft(i,1)
      do while((j.le.istop).and.(j.gt.0)) !as long as there is a child
        i=j !go one level higher
        if(i.eq.istop) exit !case of just one child
        if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
        j=ishft(i,1)
      end do

      !search correct position along the path
      do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
        i=ishft(i,-1)
      end do

      !move data
      t1=Iarray(istart)
      do while (i.ge.istart)
        t2=Iarray(i)
        Iarray(i)=t1
        t1=t2
        i=ishft(i,-1)
      end do
    end subroutine reheap

    !************************************************************************

    subroutine heapsort(Iarray)
      integer(i64),dimension(:)::Iarray
      integer(i64)::t
      integer :: i,n
      n=ubound(Iarray,1)
      ! heap creation phase (Maxheap)
      do i=ishft(n,-1),1,-1
        call reheap(Iarray,i,n)
      end do
      ! selection phase
      do i=n,2,-1
        t=Iarray(i)
        Iarray(i)=Iarray(1)
        Iarray(1)=t
        call reheap(Iarray,1,i-1)
      end do
    end subroutine heapsort

    !************************************************************************

    recursive subroutine quicksort(Iarray)
      integer(i64),dimension(:)::Iarray
      integer(i64)::t,temp
      integer::l,u,i,j
      real:: r
      integer,parameter::cutoff=50
      l=1
      u=ubound(Iarray,1)
      do while ((u-l)>=cutoff)
        ! 1.) choice of pivot
        call random_number(r)
        i=l+floor(r*(u-l+1))
        t=Iarray(i)
        Iarray(i)=Iarray(l)
        Iarray(l)=t
        ! 2.) partitioning and positioning of the pivot to position j
        i=l
        j=u
        do
          i=i+1
          do while ((i.lt.u).and.(Iarray(i).lt.t))
            i=i+1
          end do
          do while (Iarray(j).gt.t)
            j=j-1
          end do
          if (i.gt.j) exit
          temp=Iarray(i)
          Iarray(i)=Iarray(j)
          Iarray(j)=temp
          j=j-1
        end do
        Iarray(l)=Iarray(j)
        Iarray(j)=t
        ! 3.) recursion (intelligent)
        if((j-l).gt.(u-j)) then
           call quicksort(Iarray(j+1:u))
           u=j-1
        else
           call quicksort(Iarray(l:j-1))
           l=j+1
        endif
      end do
    end subroutine quicksort

    !************************************************************************

    subroutine insertsort(Iarray)
      integer(i64),dimension(:)::Iarray
      integer(i64)::t
      integer::i,j
      do i=2,ubound(Iarray,1)
        t=Iarray(i)
        j=i-1
        do while (Iarray(j)>t)
          j=j-1
	  if (j .eq. 0) exit
        end do
        Iarray(j+2:i)=Iarray(j+1:i-1)
        Iarray(j+1)=t
      end do
    end subroutine insertsort
  end subroutine sys_i64_sort


!************************************************************************


!<subroutine>
  subroutine sys_sp_sort(Darray, csortMethod)

    !<description>
    ! sort routine for single precision
    !</description>

    !<input>

    !sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer,optional::csortMethod

    !</input>

    !<inoutput>

    ! singe precision array to be sorted
    real(sp),dimension(:)::Darray
    !</inoutput>
!</subroutine>

    !if the optional argument csortMethod is present
    if(present(csortMethod)) then
       select case (csortMethod)
       case(SORT_HEAP)
          call heapsort(Darray)
       case(SORT_QUICK)
          call random_seed
          call quicksort(Darray)
          call insertsort(Darray)
       case(SORT_INSERT)
          call insertsort(Darray)
       case default
          stop 'sys_sp_sort: unknown Method'
       end select
    else
       call heapsort(Darray)
    endif
  contains

    !************************************************************************

    subroutine reheap(Darray, istart, istop)
      integer::istart,istop
      real(sp),dimension(1:istop)::Darray
      real(sp)::t1,t2
      integer::i,j
      !if(istop.eq.istart) return ! nothing to correct

      !trace the path of the bigger children
      i=istart
      j=ishft(i,1)
      do while((j.le.istop).and.(j.gt.0)) !as long as there is a child
        i=j !go one level higher
        if(i.eq.istop) exit !case of just one child
        if(Darray(i+1).gt.Darray(i)) i=i+1 !path of bigger child
        j=ishft(i,1)
      end do

      !search correct position along the path
      do while ((i.gt.istart).and.(Darray(i).lt.Darray(istart)))
        i=ishft(i,-1)
      end do

      !move data
      t1=Darray(istart)
      do while (i.ge.istart)
        t2=Darray(i)
        Darray(i)=t1
        t1=t2
        i=ishft(i,-1)
      end do
    end subroutine reheap

    !************************************************************************

    subroutine heapsort(Darray)
      real(sp),dimension(:)::Darray
      real(sp)::t
      integer :: i,n
      n=ubound(Darray,1)
      ! heap creation phase (maxheap)
      do i=ishft(n,-1),1,-1
        call reheap(Darray,i,n)
      end do
      ! selection phase
      do i=n,2,-1
        t=Darray(i)
        Darray(i)=Darray(1)
        Darray(1)=t
        call reheap(Darray,1,i-1)
      end do
    end subroutine heapsort

    !************************************************************************

    recursive subroutine quicksort(Darray)
      real(sp),dimension(:)::Darray
      real(sp)::t,temp
      integer::l,u,i,j
      real:: r
      integer,parameter::cutoff=50
      l=1
      u=ubound(Darray,1)
      do while ((u-l)>=cutoff)
        ! 1.) choice of pivot
        call random_number(r)
        i=l+floor(r*(u-l+1))
        t=Darray(i)
        Darray(i)=Darray(l)
        Darray(l)=t
        ! 2.) partitioning and positioning of the pivot to position j
        i=l
        j=u
        do
          i=i+1
          do while ((i.lt.u).and.(Darray(i).lt.t))
            i=i+1
          end do
          do while (Darray(j).gt.t)
            j=j-1
          end do
          if (i.gt.j) exit
          temp=Darray(i)
          Darray(i)=Darray(j)
          Darray(j)=temp
          j=j-1
        end do
        Darray(l)=Darray(j)
        Darray(j)=t
        ! 3.) recursion (intelligent)
        if((j-l).gt.(u-j)) then
           call quicksort(Darray(j+1:u))
           u=j-1
        else
           call quicksort(Darray(l:j-1))
           l=j+1
        endif
      end do
    end subroutine quicksort

    !************************************************************************

    subroutine insertsort(Darray)
      real(sp),dimension(:)::Darray
      real(sp)::t
      integer::i,j
      do i=2,ubound(Darray,1)
        t=Darray(i)
        j=i-1
        do while (Darray(j)>t)
          j=j-1
	  if (j .eq. 0) exit
        end do
        Darray(j+2:i)=Darray(j+1:i-1)
        Darray(j+1)=t
      end do
    end subroutine insertsort
  end subroutine sys_sp_sort


!************************************************************************


!<subroutine>
  subroutine sys_dp_sort(Darray, csortMethod, Imapping)

    !<description>
    !Sorting routine for double precision arrays. If more than one vector must be
    !sorted, or if several vectors may be sorted the same way, the mapping
    !array Imapping can be computed. This array must be created outside the routine
    !and must have the same length as Iarray. Then, the other vectors are computed
    !by Vec(Imapping(i)). The array Imapping is an optional argument. The sorting
    !method used is told by the parameter csortMethod. If this optional argument
    !is not used, the routine performs heapsort.
    !</description>

    !<input>

    !sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer,optional::csortMethod

    !</input>

    !<inoutput>

    ! double precision array to be sorted
    real(dp),dimension(:)::Darray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer,dimension(:), optional:: Imapping

    !</inoutput>
!</subroutine>

    !if the optional argument csortMethod is present
    if(present(csortMethod)) then
      select case (csortMethod)
      case(SORT_HEAP)
        if (present(Imapping)) then
          call heapsort(Darray, Imapping)
        else
          call heapsort(Darray)
        endif
      case(SORT_QUICK)
        call random_seed
        
        if (present(Imapping)) then
          call quicksort(Darray, Imapping)
          call insertsort(Darray, Imapping)
        else
          call quicksort(Darray)
          call insertsort(Darray)
        endif
      case(SORT_INSERT)
        
        if (present(Imapping)) then
          call insertsort(Darray, Imapping)
        else
          call insertsort(Darray)
        endif
        
      case default
        stop 'sys_dp_sort: unknown Method'
      end select
    else
      
      if (present(Imapping)) then
        call heapsort(Darray, Imapping)
      else
        call heapsort(Darray)
      endif
      
    endif
    
  contains

    !************************************************************************

    subroutine reheap(Darray, istart, istop, Imapping)
      integer::istart,istop
      real(dp),dimension(1:istop)::Darray
      integer,dimension(1:istop), optional:: Imapping
      real(dp)::t1,t2
      integer :: i1, i2
      integer::i,j
      !if(istop.eq.istart) return ! nothing to correct

      if (present(Imapping)) then

        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j.gt.0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Darray(i+1).gt.Darray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        end do
        
        !search correct position along the path
        do while ((i.gt.istart).and.(Darray(i).lt.Darray(istart)))
          i=ishft(i,-1)
        end do
        
        !move data
        t1=Darray(istart)
        i1=Imapping(istart)

        do while (i.ge.istart)
          t2=Darray(i)
          Darray(i)=t1
          t1=t2

          i2=Imapping(i)
          Imapping(i)=i1
          i1=i2
          i=ishft(i,-1)
        end do

      else

        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j.gt.0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Darray(i+1).gt.Darray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        end do
        
        !search correct position along the path
        do while ((i.gt.istart).and.(Darray(i).lt.Darray(istart)))
          i=ishft(i,-1)
        end do
        
        !move data
        t1=Darray(istart)
        do while (i.ge.istart)
          t2=Darray(i)
          Darray(i)=t1
          t1=t2
          i=ishft(i,-1)
        end do
      endif
    end subroutine reheap

    !************************************************************************

    subroutine heapsort(Darray, Imapping)
      real(dp),dimension(:)::Darray
      integer,dimension(:), optional:: Imapping
      real(dp)::t
      integer :: t2
      integer :: i,n
      n=ubound(Darray,1)
      ! heap creation phase (maxheap)

      if (present(Imapping)) then
        do i=ishft(n,-1),1,-1
          call reheap(Darray,i,n, Imapping)
        end do
        ! selection phase
        do i=n,2,-1
          t2=Imapping(i)
          Imapping(i)=Imapping(1)
          Imapping(1)=t2

          t=Darray(i)
          Darray(i)=Darray(1)
          Darray(1)=t
          call reheap(Darray,1,i-1, Imapping)
        end do
      else

        do i=ishft(n,-1),1,-1
          call reheap(Darray,i,n)
        end do
        ! selection phase
        do i=n,2,-1
          t=Darray(i)
          Darray(i)=Darray(1)
          Darray(1)=t
          call reheap(Darray,1,i-1)
        end do
      endif
    end subroutine heapsort

    !************************************************************************

    recursive subroutine quicksort(Darray, Imapping)
      real(dp),dimension(:)::Darray
      integer,dimension(:), optional:: Imapping
      real(dp)::t, temp
      integer :: t2, temp2
      integer::l,u,i,j
      real:: r
      integer,parameter::cutoff=50
      l=1
      u=ubound(Darray,1)
      
      if (present(Imapping)) then
        do while ((u-l)>=cutoff)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Darray(i)
          Darray(i)=Darray(l)
          Darray(l)=t

          t2=Imapping(i)
          Imapping(i)=Imapping(l)
          Imapping(l)=t2

          ! 2.) partitioning and positioning of the pivot to position j
          do
            i=i+1
            do while ((i.lt.u).and.(Darray(i).lt.t))
              i=i+1
            end do
            do while (Darray(j).gt.t)
              j=j-1
            end do
            if (i.gt.j) exit
            temp=Darray(i)
            temp2=Imapping(i)

            Darray(i)=Darray(j)
            Darray(j)=temp

            Imapping(i)=Imapping(j)
            Imapping(j)=temp2

            j=j-1
          end do

          Darray(l)=Darray(j)
          Darray(j)=t
          
          Imapping(l)=Imapping(j)
          Imapping(j)=t2

          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Darray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            call quicksort(Darray(l:j-1), Imapping(l:j-1))
            l=j+1
          endif
        end do
        
      else
        do while ((u-l)>=cutoff)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Darray(i)
          Darray(i)=Darray(l)
          Darray(l)=t
          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Darray(i).lt.t))
              i=i+1
            end do
            do while (Darray(j)>t)
              j=j-1
            end do
            if (i.gt.j) exit
            temp=Darray(i)
            Darray(i)=Darray(j)
            Darray(j)=temp
            j = j-1
          end do

          Darray(l)=Darray(j)
          Darray(j)=t
          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Darray(j+1:u))
            u=j-1
          else
            call quicksort(Darray(l:j-1))
            l=j+1
          endif
        end do
      endif

    end subroutine quicksort

    !************************************************************************

    subroutine insertsort(Darray, Imapping)
      real(dp),dimension(:)::Darray
      integer,dimension(:), optional:: Imapping

      real(dp)::t
      integer :: t2
      integer::i,j

      if (present(Imapping)) then
        do i=2, ubound(Darray,1)
          t=Darray(i)
          t2 = Imapping(i)
          j=i-1
          do while (Darray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          end do
          Darray(j+2:i)=Darray(j+1:i-1)
          Darray(j+1)=t

          Imapping(j+2:i)=Imapping(j+1:i-1)
          Imapping(j+1)=t2

        end do
      else
        do i=2,ubound(Darray,1)
          t=Darray(i)
          j=i-1
          do while (Darray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          end do
          Darray(j+2:i)=Darray(j+1:i-1)
          Darray(j+1)=t
      end do
    endif
  end subroutine insertsort
end subroutine sys_dp_sort


!*****************************************************************************************


!<function>
  real(DP) function sys_triArea(dx1,dy1,dx2,dy2,dx3,dy3)

    !<description>
    !Compute the oriented area of the triangle spanned by (dxi, dyi).
    !</description>

    !<result>
    !area of the triangle
    !</result>

    !<input>

    !coordinates of the three points
    real(DP) :: dx1,dy1,dx2,dy2,dx3,dy3

    !</input>
!</function>

    sys_triArea = 0.5*( (dx1-dx3)*(dy2-dy3) - (dy1-dy3)*(dx2-dx3))

  end function sys_triArea


!*****************************************************************************************


!<function>
  real(DP) function sys_quadArea(dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4)

    !<description>
    !Compute the oriented area of the quadrilateral spanned by (dxi, dyi).
    !</description>

    !<result>
    !area of the quadrilateral
    !</result>

    !<input>

    !coordinates of the four points
    real(DP) :: dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4

    !</input>
!</function>

    sys_quadArea = 0.5*( (dx1*dy2 - dy1*dx2) + (dx2*dy3 - dy2*dx3) &
                         +(dx3*dy4 - dy3*dx4) + (dx4*dy1 - dy4*dx1))

  end function sys_quadArea


!*****************************************************************************************


!<subroutine>
    subroutine sys_parInElement(DcoordX, DcoordY, dxpar, dypar, dxreal, dyreal)

    !<description>
    !This subroutine is to find the parameter values for a given point (x,y) in real
    !coordinates.\\
    !
    !Remark: This is a difficult task, as usually in FEM codes the parameter values are
    !known and one wants to obtain the real coordinates.

    !inverting the bilinear trafo in a straightforward manner by using pq-formula
    !does not work very well, as it is numerically unstable. For parallelogram-shaped
    !elements, one would have to introduce a special treatment. For nearly
    !parallelogram-shaped elements, this can cause a crash as the argument of the
    !square root can become negative due to rounding errors. In the case of points near
    !the element borders, we divide nearly 0/0.\\
    !
    !Therefore, we have implemented the algorithm described in
    ! Introduction to Fineite Element Methods, Carlos Felippa, Department of Aerospace
    ! Engineering Sciences and Center for Aerospace Structures,
    ! http://titan.colorado.edu/courses.d/IFEM.d/  .
    !</description>

!<input>

    !coordinates of the evaluation point
    real(DP),intent(in) :: dxreal, dyreal

    !coordinates of the element vertices
    real(DP), dimension(4),intent(in) :: DcoordX
    real(DP), dimension(4),intent(in) :: DcoordY

!</input>

!<output>

    !parameter values of (x,y)
    real(DP), intent(out) :: dxpar,dypar

!</output>

!</subroutine>

    !internal variables
    real(DP) :: x1,x2,x3,x4,y1,y2,y3,y4,xb,yb,xcx,ycx,xce,yce
    real(DP) :: A,J1,J2,x0,y0,bxi,beta,cxi,xp0,yp0,ceta,root1,root2
    real(DP) :: xip1,xip2,etap1,etap2,d1,d2


    !Get nodal x-coordinates
    X1 = DcoordX(1)
    X2 = DcoordX(2)
    X3 = DcoordX(3)
    X4 = DcoordX(4)

    !Get nodal y-coordinates
    Y1 = DcoordY(1)
    Y2 = DcoordY(2)
    Y3 = DcoordY(3)
    Y4 = DcoordY(4)

    XB = X1-X2+X3-X4
    YB = Y1-Y2+Y3-Y4

    XCX = X1+X2-X3-X4
    YCX = Y1+Y2-Y3-Y4

    XCE = X1-X2-X3+X4
    YCE = Y1-Y2-Y3+Y4

    A = 0.5*((X3-X1)*(Y4-Y2)-(X4-X2)*(Y3-Y1))

    J1 = (X3-X4)*(Y1-Y2)-(X1-X2)*(Y3-Y4)
    J2 = (X2-X3)*(Y1-Y4)-(X1-X4)*(Y2-Y3)

    X0 = 0.25*(X1+X2+X3+X4)
    Y0 = 0.25*(Y1+Y2+Y3+Y4)

    XP0 = dxreal-X0
    YP0 = dyreal-Y0

    BXI  =  A-XP0*YB+YP0*XB
    BETA = -A-XP0*YB+YP0*XB

    CXI  = XP0*YCX-YP0*XCX
    CETA = XP0*YCE-YP0*XCE

    ROOT1 = -SQRT(BXI**2-2.0*J1*CXI)-BXI
    ROOT2 =  SQRT(BXI**2-2.0*J1*CXI)-BXI
    IF (ROOT1.NE.0D0) THEN
      XIP1 = 2.0*CXI/ROOT1
    ELSE
      XIP1 = 1E15
    ENDIF
    IF (ROOT2.NE.0D0) THEN
      XIP2 = 2D0*CXI/ROOT2
    ELSE
      XIP2 = 1E15
    ENDIF

    ROOT1 =  SQRT(BETA**2+2D0*J2*CETA)-BETA
    ROOT2 = -SQRT(BETA**2+2D0*J2*CETA)-BETA
    IF (ROOT1.NE.0D0) THEN
      ETAP1 = 2D0*CETA/ROOT1
    ELSE
      ETAP1 = 1D15
    ENDIF
    IF (ROOT2.NE.0D0) THEN
      ETAP2 = 2D0*CETA/ROOT2
    ELSE
      ETAP2 = 1D15
    ENDIF

    D1 = SQRT(XIP1**2+ETAP1**2)
    D2 = SQRT(XIP2**2+ETAP2**2)

    IF (D1.LT.D2) THEN
      dxpar = XIP1
      dypar = ETAP1
    ELSE
      dxpar = XIP2
      dypar = ETAP2
    ENDIF

  end subroutine sys_parInElement


!************************************************************************************


!<function>
  integer function ipointInElement(DcoordX, DcoordY, dx,dy)
!<description>
!This function is to detect wether a given point (dx, dy) is inside the macro given
!by the vertices coordinates in Dcoord. We devide the macro in  4 triangles.
!The point P lies in the quadrilateral ABCD, if the oriented areas of ABP,BCP,CDP,DAP
!all are lower than 0. The area is computed in the function area. This requires the
!counterclockwise numeration of the quad.
!</description>

!<input>

    !coordinates of vertices of element
    real(DP), dimension(4) :: DcoordX
    real(DP), dimension(4) :: DcoordY

    !coordinates of test point
    real(DP) :: dx, dy
!</input>

!<result>
    !        -i, if point is vertex i\\
    !         0, if point is outside the macro\\
    !         1, if point is in macro or on a macro edge or is a vertex\\
    !         2, if point is in the interior of the macro\\
    !
!</result>
!</function>
    integer :: i                                !index variable


    ipointInElement = 0
    do i = 1,4

       if ((DcoordX(i)-dx)*(DcoordX(i)-dx) + &
            (DcoordY(i)-dy)*(DcoordY(i)-dy) .lt. 1d-14) then
          ipointInElement = -i
          GOTO 99999
       endif
    enddo

    if ((sys_triArea(DcoordX(1),DcoordY(1),DcoordX(2),DcoordY(2),dx,dy)     &
         .ge.(-1.0D-14)) .and.                                                  &
        (sys_triArea(DcoordX(2),DcoordY(2),DcoordX(3),DcoordY(3),dx,dy)     &
         .ge.(-1.0D-14)) .and.                                                  &
        (sys_triArea(DcoordX(3),DcoordY(3),DcoordX(4),DcoordY(4),dx,dy)     &
         .ge.(-1.0D-14)) .and.                                                  &
        (sys_triArea(DcoordX(4),DcoordY(4),DcoordX(1),DcoordY(1),dx,dy)     &
         .ge.(-1.0D-14))                                                        &
       ) then
       ipointInElement = 1
    endif

    if ((sys_triArea(DcoordX(1),DcoordY(1),DcoordX(2),DcoordY(2),dx,dy)     &
         .ge.(1.0D-14)).and.                                                    &
        (sys_triArea(DcoordX(2),DcoordY(2),DcoordX(3),DcoordY(3),dx,dy)     &
         .ge.(1.0D-14)).and.                                                    &
        (sys_triArea(DcoordX(3),DcoordY(3),DcoordX(4),DcoordY(4),dx,dy)     &
         .ge.(1.0D-14)).and.                                                    &
        (sys_triArea(DcoordX(4),DcoordY(4),DcoordX(1),DcoordY(1),dx,dy)     &
         .ge.(1.0D-14))                                                         &
       ) then
       ipointInElement = 2
    endif
99999 end function ipointInElement

end module

