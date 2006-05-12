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
!# <name> output </name>                                                #
!#                                                                      #
!#                                                                      #
!# <purpose>                                                            #
!# Handling the output of messages to terminals and logfiles            #
!# </purpose>                                                           #
!#                                                                      #
!########################################################################

!# Current version: $Id: output.f90,v 1.23 2006/03/22 15:12:44 goeddeke4 Exp $

!<!--
!#include <feastdefs.h>
! -->

module genoutput

  use fsystem

  implicit none

! <constants>
! <constantblock>

  ! device number for logfile
  integer :: OU_LOG

  ! device number for I/O
  integer, parameter :: OU_CHANNEL1 = 32

  ! device number for I/O
  integer, parameter :: OU_CHANNEL2 = 33

  ! device number for I/O
  integer, parameter :: OU_CHANNEL3 = 34

  integer, parameter :: OU_CHANNEL_FEAST=35
  integer, parameter :: OU_CHANNEL_FGEO=36
  integer, parameter :: OU_CHANNEL_FFGEO=37
  integer, parameter :: OU_CHANNEL_FMESH=38
  integer, parameter :: OU_CHANNEL_FPART=39
  integer, parameter :: OU_CHANNEL_FBC=40

  ! device number for result files
  integer, parameter :: OU_CHANNEL_RESULT=41

  ! output_cmode, trace info in log file
  integer, parameter :: OM_LOG = 0

  ! output_cmode, trace info on terminal
  integer, parameter :: OM_TERM = 1

  ! output_cmode, trace info in file/term
  integer, parameter :: OM_BOTH = 2

  ! output_cmode, no trace information
  integer, parameter :: OM_NONE = 3

  ! output_clevel, no trace information
  integer, parameter :: OL_NONE = 0


  ! output_clevel, only output of
  ! error messages
  integer, parameter :: OL_ERROR = 1

  ! output_clevel, only system messages
  integer, parameter :: OL_MSG = 2

  ! output_clevel, only system/timer messages
  integer, parameter :: OL_TIMER = 3

  ! output_clevel, trace info level 1
  integer, parameter :: OL_TRACE1 = 4

  ! output_clevel, trace info level 2
  integer, parameter :: OL_TRACE2 = 5

  ! output_clevel, trace info level 3
  integer, parameter :: OL_TRACE3 = 6
! </constantblock>
! </constants>


! <globals>

  ! level of tracing information
  integer :: output_clevel

  ! number of actual result file
  integer :: output_numberResultFile
  
! </globals>


  integer :: ioutputoffset

  contains


!************************************************************************

!<function>
  integer function output_cOutputLevel()

!<description>
! This function returns the output / verbosity level.
!</description>

! </function>
    
     output_cOutputLevel = output_clevel
  
  end function output_cOutputLevel


!************************************************************************

!<subroutine>
  subroutine output_openResultFile()

!<description>
!This subroutine opens the result file for global result data. 
!</description>

!</subroutine>

    character (len=SYS_STRLEN) :: sfilename
    character (len=SYS_STRLEN) :: stmp
    integer :: iostat
    
    write (stmp, '(I16)') output_numberResultFile
    
    sfilename = "result.txt"
    
    open(unit = OU_CHANNEL_RESULT, file = trim(sfilename), IOSTAT = iostat, ACTION = 'WRITE')
    
    if (iostat .ne. 0) then
      OU_LOG        = 6
      output_clevel = OL_ERROR
      call output_line(OL_ERROR, "", "Unable to open FEAST result file named " // &
                                     trim(sfilename) // ",")
      call output_line(OL_ERROR, "", "aborting program.")
      write(OU_LOG, *) "FATAL EXIT!"
      call output_close
      !The problem is that in this module nothing is known about
      !a parallel environment and thus, the parallel process hangs now
      !indefinitely!
      stop       
    endif
    
  end subroutine output_openResultFile
  

!************************************************************************
  
!<subroutine>
  subroutine output_closeResultFile()

!<description>
!This routine closes the result file.
!</description>

!</subroutine>
 
    call sys_flush(OU_CHANNEL_RESULT)

    open(unit = OU_CHANNEL_RESULT)
    
    output_numberResultFile = output_numberResultFile + 1
      
  end subroutine output_closeResultFile

  
!************************************************************************

!<subroutine>
  subroutine output_writeIntResult(stag, ivalue)
  
!<description>
!This routine writes an integer formatted result item to the result file
!</description>

!<input>
    
    !tag name of the result data 
    character(len=*) :: stag
    
    !result data
    integer :: ivalue

!</input>

!</subroutine>
    character (len=10) :: sstring

    sstring = sys_siL(ivalue, 10)
    write (OU_CHANNEL_RESULT, "(A,A,A)") stag, " = ", trim(sstring)
     
  end subroutine output_writeIntResult


!************************************************************************

!<subroutine>
  subroutine output_writeDoubleResult(stag, dvalue)

!<description>
!This routine writes an double formatted result item to the result file
!</description>

!<input>

    !tag name of the result data 
    character(len=*) :: stag

    !result data
    real(DP) :: dvalue

!</input>

!</subroutine>
    character (len=16) :: sstring

    sstring = trim(sys_sdEL(dvalue, 8))
    write (OU_CHANNEL_RESULT, "(A,A,A)") stag, " = ", trim(sstring)
    
  end subroutine output_writeDoubleResult


!************************************************************************

!<subroutine>
  subroutine output_writeFloatResult(stag, dvalue)

!<description>
!This routine writes an single formatted result item to the result file
!</description>

!<input>

    !tag name of the result data 
    character(len=*) :: stag

    !result data
    real(DP) :: dvalue

!</input>

!</subroutine>

    character (len=10) :: sstring

    sstring = sys_sdL(dvalue, 3)
    write (OU_CHANNEL_RESULT, "(A,A,A)") stag, " = ", trim(sstring)

  end subroutine output_writeFloatResult


!************************************************************************

!<subroutine>
  subroutine output_writeStringResult(sstring)

!<description>
!This routine writes a string to the result file
!</description>

!<input>
    !result data 
    character(len=*) :: sstring
!</input>

!</subroutine>

    write (OU_CHANNEL_RESULT, "(A)") trim(sstring)

  end subroutine output_writeStringResult
  
     
!************************************************************************

!<function>
  logical function output_do(cpriorityLevel)

!<description>
! This function returns the do state of a pending output operation
!</description>

! <input>
    !priority level, OL_
    integer :: cpriorityLevel
! </input>

! </function>
    
    output_do=  (output_clevel .ge. cpriorityLevel)
  
  end function output_do


!************************************************************************

!<subroutine>
  subroutine output_line(cpriorityLevel, sroutine, smessage)

!<description>
! This routine writes a message with priority cpriorityLevel to log device
!</description>

! <input>
    !priority level, OL_
    integer :: cpriorityLevel

    !name of calling routine
    character (len=*) :: sroutine

    !message which shall be printed
    character (len=*) :: smessage

! </input>

! </subroutine>

    character (len=64) :: sbar


    if (output_clevel .ge. cpriorityLevel) then

      select case (cpriorityLevel)
        case (OL_TRACE1)
          sbar = "*** " // sroutine // ":"
        case (OL_TRACE2)
          sbar = "****** " // sroutine//":"
        case (OL_TRACE3)
          sbar = "********* " // sroutine // ":"
        case default
          sbar = "* " // sroutine // ":"
      end select

      if (OU_LOG .eq. 6) then
        write(6,'(A,A,A)') trim(sbar), " ", trim(smessage)
      else
        write(OU_LOG, '(A,A,A)') trim(sbar), " ", trim(smessage)
        if (cpriorityLevel .eq. OL_ERROR) then
          write(6,'(A,A,A)') trim(sbar), " ", trim(smessage)
        endif
      endif
    endif

  end subroutine output_line

!************************************************************************

!<subroutine>
  subroutine output_init(shome, cpriorityLevel)

!<description>
! This subroutine initializes the log files.
!</description>

! <input>
    !priority level
    integer :: cpriorityLevel

    ! path to home directory
    character (len=*) :: shome
! </input>

! </subroutine>

    !complete name of the FEAST-log file (e.g. FEASTLOG.001)
    character (len=256) :: sfilename

    !index of calling parallel block (this time as character)
    character (len=3) :: stmp

    !return value for 'file open' statement
    integer :: iostat = 0

    ! Initialise the output device
    OU_LOG      = 48

    sfilename = trim(shome) // "/feastlog." // trim(stmp)

    open(unit = OU_LOG, file = trim(sfilename), IOSTAT = iostat, ACTION = 'WRITE')

    output_numberResultFile=0
    
    if (iostat .ne. 0) then
      OU_LOG        = 6
      output_clevel = OL_ERROR
      call output_line(OL_ERROR, "", "Unable to open FEAST log file named " // &
                                     trim(sfilename) // ",")
      call output_line(OL_ERROR, "", "aborting program.")
      write(OU_LOG, *) "FATAL EXIT!"
      call output_close
      !The problem is that in this module nothing is known about
      !a parallel environment and thus, the parallel process hangs now
      !indefinitely!
      stop
    endif

    output_clevel = cpriorityLevel

  end subroutine output_init

!************************************************************************

!<subroutine>
  subroutine output_close()

!<description>
! This routine closes the log device.
!</description>

! </subroutine>

    close(OU_LOG)
    
  end subroutine output_close


end module genoutput
!replacements to be done:
!   44 OM_LOG                                  ->  OU_MODE_LOG  (etc.)
!   49 OL_NONE                                 ->  OU_LEV_NONE  (etc.)
!  160 cmodus                                  ->  ctimerModus
!  161 rt                                      ->  rtimeStamp
