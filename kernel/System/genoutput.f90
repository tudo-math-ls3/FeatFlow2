!##############################################################################
!# ****************************************************************************
!# <name> genoutput </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains several routines for input/output of messages to the
!# terminal and/or to log files.
!#
!# The following routines can be found here:
!#
!# 1.) output_init
!#     -> Initialises the output. Opens a log file and assigns it an
!#        output channel OU_LOG.
!#
!# 2.) output_line
!#     -> Writes a message to the terminal and/or the log file.
!#
!# 3.) output_separator
!#     -> Writes a separator line to the terminal and/or the log file.
!#
!# 4.) output_lbrk
!#     -> Writes an empty line / line break to the terminal and/or the log 
!#        file.
!#
!# 5.) output_simple
!#     -> Writes a message to the terminal and/or the log file.
!#        For old FEAT compatibility. Uses a priority identifier to decide on
!#        whether to write only to the log file or to both, log file
!#        and terminal.
!#
!# 6.) output_simple_sep
!#     -> Writes a separator line to the terminal and/or the log file.
!#        For old FEAT compatibility. Uses a priority identifier to decide on
!#        whether to write only to the log file or to both, log file
!#        and terminal.
!#
!# 7.) output_done
!#     -> Closes the log file, releases all ressources in use.
!#
!# HOW TO USE:
!# -----------
!# As long as output_init is not called, all output is redirected only to the 
!# standard output channel. Directly after program start, output_init should
!# be called so that the log file is written.
!#
!# Proceed as follows:
!#
!# 1.) CALL output_init ()
!#     -> Initialise the output system for output on the terminal
!#
!#     Alternatively one can use:
!#
!#     CALL output_init ('mylogfile.txt')
!#     -> Opens a log file 'mylogfile.txt' for the output additional to
!#        terminal output
!#
!#     Alternatively one can use
!# 
!#     CALL output_init ('mylogfile.txt','myerrorlogfile.txt')
!#     -> Opens a log file 'mylogfile.txt' for the output and
!#        'myerrorlogfile.txt' for error output (both additionally to
!#        terminal output).
!#     
!# 2.) CALL output_line ('This is a message')
!#
!#     -> Writes a message to the terminal. If output_init was used
!#        before, output_line_std will also write the message to the log file.
!#
!#     Alternatively, one can specify where to write output to by using
!#     different variants of output_line:
!#
!#     a) CALL output_lbrk ()
!#
!#     -> Writes an empty line to the terminal and log file. This is the same
!#        as CALL output_line ('').
!#
!#     b) CALL output_line('A message only to the terminal.', &
!#                         OU_CLASS_MSG,OU_MODE_TERM)
!#
!#     -> Writes a message only to the terminal
!#
!#     c) CALL output_line ('A message only to the log file.', &
!#                          OU_CLASS_MSG,OU_MODE_LOG)
!#
!#     -> Writes a message only to the log file
!#
!#     d) CALL output_line ('A special debug message.', &
!#                          OU_CLASS_TRACE1,OU_MODE_STD,'mysubroutine')
!#
!#     -> Writes a debug message with '*** (mysubroutine):' in front to the 
!#        terminal and the log file. This is usually used for debug purposes.
!#
!#     e) CALL output_line ('This is an error message.', &
!#                          OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
!#
!#     -> Writes a debug message with 'Error (mysubroutine):' in front to the 
!#        error log file and error output channel.
!#
!#       or even simpler:
!#
!#        CALL output_line ('This is an error message.',OU_CLASS_ERROR)
!#        CALL output_line ('This is an warning message.',OU_CLASS_WARNING)
!#
!#     f) CALL output_separator (OU_SEP_MINUS)
!#
!#     -> Writes a separation line with '-' signs to the terminal / log file
!#
!#     g) MT = 1
!#        CALL output_simple (MT,'A log file message.')
!#     
!#     -> Writes a message to the log file, not to the terminal. FEAT1.0
!#        compatibility routine for "MT=1"
!#
!#     h) MT = 2
!#        CALL output_simple (MT,'A log file message')
!#     
!#     -> Writes a message to the terminal and to the log file. FEAT1.0
!#        compatibility routine for "MT=2".
!#
!#     i) MT = 2
!#        CALL output_simple (MT)
!#     
!#     -> Writes an empty line to the terminal and to the log file. FEAT1.0
!#        compatibility routine for "MT=2".
!#
!#     i) MT = 2
!#        CALL output_simple_sep (MT,OU_SEP_MINUS)
!#     
!#     -> Writes a separation line with '-' signs to the terminal and 
!#        to the log file. FEAT1.0 compatibility routine for "MT=2".
!#
!#     For futher possibilities, consider the documentation out output_line.
!#
!# 3.) CALL output_done()
!#
!#     -> Closes the output channel(s).
!#
!# Logging benchmark data
!# ----------------------
!# The output library furthermore supports the (semi-)automatic output of
!# deterministic data to a so callen 'benchmark log file'. This is typically
!# used to write out data to a specific file which is compared to reference
!# data in regression tests. Such data is meant to be deterministic,
!# reproducable in every run and usually formatted in a 'nice' style that
!# differences to reference results can easily be checked.
!#
!# To activate the benchmark log file, one has to specify an additional
!# parameter in the call to output_init:
!#
!#    CALL output_init ('mylogfile.txt','myerrorlogfile.txt','benchmarkresultfile')
!#
!# This opens a file 'benchmarkresultfile' where benchmark data is written to.
!# To write a string to this file, the application has to call output_line
!# with an extended output-mode:
!#
!#   CALL output_line('A message only to the benchmark log.', &
!#                    OU_CLASS_MSG,OU_MODE_BENCHLOG)
!#
!# This writes a message directly to the benchmark log file. It's also possible
!# to write out data to the standard terminal plus the benchmark log file
!# by combining the output mode constants:
!#
!#   CALL output_line('A message only to the benchmark log.', &
!#                    OU_CLASS_MSG,OU_MODE_STD+OU_MODE_BENCHLOG)
!#
!# In all cases, the additional constant OU_MODE_BENCHLOG must be manually
!# specified as logging to the benchmark log file is not done automatically.
!#
!# </purpose>
!##############################################################################

module genoutput

  use fsystem

  implicit none

!<constants>

!<constantblock description="Output mode. Decides on where the output is written to.">

  ! Output mode: Write to main log file / error log file (depending on whether
  ! a message is considered as an error or not by OU_CLASS_ERROR/OU_CLASS_WARNING)
  integer(I32), parameter :: OU_MODE_LOG      = 2**0
  
  ! Output mode: Write to terminal
  integer(I32), parameter :: OU_MODE_TERM     = 2**1
  
  ! Output mode: Write to benchmark log file
  integer(I32), parameter :: OU_MODE_BENCHLOG = 2**2

  ! Output mode: Write to both, log file and terminal
  integer(I32), parameter :: OU_MODE_STD      = OU_MODE_LOG+OU_MODE_TERM
  
!</constantblock>
  
!<constantblock description="Output classification. Prints an additional classification string\\
!                            when writing a string to a log file / terminal">

  ! Output classification: Standard message
  integer, parameter :: OU_CLASS_MSG      = 0
  
  ! Output classification: Trace information, level 1
  integer, parameter :: OU_CLASS_TRACE1   = 1

  ! Output classification: Trace information, level 2
  integer, parameter :: OU_CLASS_TRACE2   = 2

  ! Output classification: Trace information, level 3
  integer, parameter :: OU_CLASS_TRACE3   = 3

  ! Output classification: System/Timer message
  integer, parameter :: OU_CLASS_SYSTEM   = 4
  
  ! Output classification: Error message.
  integer, parameter :: OU_CLASS_ERROR    = 5

  ! Output classification: Warning message.
  integer, parameter :: OU_CLASS_WARNING  = 6

!</constantblock>

!<constantblock description="Type of separator line">

  ! Separator line: MINUS character
  integer, parameter :: OU_SEP_MINUS  = 0

  ! Separator line: STAR character 
  integer, parameter :: OU_SEP_STAR   = 1

  ! Separator line: EQUAL character
  integer, parameter :: OU_SEP_EQUAL  = 2

  ! Separator line: DOLLAR character
  integer, parameter :: OU_SEP_DOLLAR = 3

  ! Separator line: @ character
  integer, parameter :: OU_SEP_AT     = 4

  ! Separator line: + character
  integer, parameter :: OU_SEP_PLUS   = 5

  ! Separator line: ~ character
  integer, parameter :: OU_SEP_TILDE  = 6

  ! Separator line: & character
  integer, parameter :: OU_SEP_AMPAND = 7

  ! Separator line: % character
  integer, parameter :: OU_SEP_PERC   = 8

  ! Separator line: # character
  integer, parameter :: OU_SEP_HASH   = 9
  
!</constantblock>

!<constantblock>

  ! Length of a line on the terminal / in the log file.
  ! Standard value = 80 characters
  integer :: OU_LINE_LENGTH         = 80

  ! Global device number for terminal output
  integer :: OU_TERMINAL            = 6

  ! Global device number for terminal input
  integer :: IN_TERMINAL            = 5

  ! Global device number for errors on terminal
  integer :: OU_ERROR               = 6
  
  ! Global device number for log file. The log file is opened in output_init.
  ! <=0: no log file output.
  integer :: OU_LOG                 = 0
  
  ! Global device number for error log file. The error log file is opened 
  ! in output_init and usually coincides with OU_LOG.
  ! <=0: write error messages to standard log file.
  integer :: OU_ERRORLOG            = 0

  ! Global device number for benchmark log file. The error log file is opened 
  ! in output_init and usually coincides with OU_LOG.
  ! <=0: write error messages to standard log file.
  integer :: OU_BENCHLOG            = 0

!</constantblock>

!</constants>  

  interface output_line
    module procedure output_line_std
    module procedure output_line_feast
  end interface
  
  interface output_init
    module procedure output_init_simple
    module procedure output_init_logfile
    module procedure output_init_standard
  end interface

contains

!************************************************************************

!<subroutine>
  subroutine output_openLogfile(sfilename, iunit)

  !<description>
    ! This routine opens a file sfilename for writing. This is 
    ! usually used to create a log file. The routine will return a
    ! handle iunit with an output channel where text can be written to.
    ! If the file exists, it's overwritten.
    !
    ! As long as this routine is not called, all output is redirected
    ! only to the standard output channel.
  !</description>

  !<input>

    ! Name of the file to overwrite.
    character(*), intent(IN) :: sfilename

  !</input>

  !<output>
    !unit of the opened file
    integer, intent(OUT) :: iunit
  !</output>
  
!</subroutine>

    integer :: istatus ! status variable for opening procedure

    if (trim(sfilename) .eq. "") then 
      write (*,'(A)') 'Error: io_openFileForWriting. sfilename undefined!'
      return
    endif

    iunit = sys_getFreeUnit()
    if (iunit .eq. -1) then
      write (*,'(A)') 'Error: io_openFileForWriting. No output channel available!'
      return
    endif

    open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
          action="write")

    if (istatus .ne. 0) then
      write (unit=*, fmt=*) 'Error: io_openFileForWriting. &
                            &Error while opening file "',&
                            trim(sfilename), '". ***'
      iunit = -1
    endif

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_init_simple ()

!<description>
  ! Initialises the output system. The output system is configured to show all
  ! messages on the terminal, no log file is opened.
!</description>
  
!</subroutine>

    call output_init_standard ("","")

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_init_logfile (slogFilename)

!<description>
  ! Initialises the output system. 
  ! If sfilename is given, a log file with that name is opened and all log and error
  ! messages are written to it. If not given, all output is directed to the 
  ! terminal only.
!</description>
  
!<input>

  ! Name of a log file for standard messages. If "" is specified, 
  ! output messages are written to the standard output.
  character(LEN=*), intent(IN) :: slogFilename
  
!</input>

!</subroutine>

    call output_init_standard (slogFilename,"")

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_init_standard (slogFilename,serrorFilename,sbenchLogfile)

!<description>
  ! Initialises the output system. 
  ! If sfilename is given, a log file with that name is opened and all log 
  ! messages are written to it. If not given, all output is directed to the 
  ! terminal only.
  ! If serrorFilename is given, an error log file with filename serrorFilename
  ! is opened and all error messages are redirected to it. If not given,
  ! the error messages are directed to the log file -- or to the terminal, if
  ! slogFilename does not exist.
!</description>
  
!<input>

  ! Name of a log file for standard messages. If ""
  ! is specified, output messages are written to the standard output.
  character(LEN=*), intent(IN) :: slogFilename
  
  ! Name of an log file for error messages. If ""
  ! is specified, errors are written to the standard log file. The name of 
  ! the file may also coincide with slogFilename.
  character(LEN=*), intent(IN) :: serrorFilename

  ! OPTIONAL: Name of an log file for deterministic benchmark messages. If 
  ! not present or set to "", benchmark messages are not written out. The 
  ! name of the file may also coincide with slogFilename.
  character(LEN=*), intent(IN), OPTIONAL :: sbenchLogfile
  
!</input>

!</subroutine>

    ! Close previously opened log files.
    call output_done ()    
    
    ! Name of the standard logfile given?
    if (slogFilename .ne. "") then
    
      ! Open a log file
      call output_openLogfile(slogFilename, OU_LOG)
    
    end if

    ! Name of the error logfile given?
    if (serrorFilename .ne. "") then
    
      ! Both filenames the same?
      if (slogFilename .eq. serrorFilename) then
        OU_ERRORLOG = OU_LOG
      else
        ! Open an error log file
        call output_openLogfile(serrorFilename, OU_ERRORLOG)
      end if
      
    else
    
      ! Write error messages to standard log file.
      OU_ERRORLOG = OU_LOG
    
    end if

    ! Name of the benchmark logfile given?
    if (present(sbenchLogfile)) then
      if (sbenchLogfile .ne. "") then
      
        ! Both filenames the same?
        if (slogFilename .eq. serrorFilename) then
          OU_BENCHLOG = OU_LOG
        else
          ! Open a benchmark log file
          call output_openLogfile(sbenchLogfile, OU_BENCHLOG)
        end if
        
      else
      
        ! No benchmark output
        OU_BENCHLOG = 0
      
      end if
      
    else

      ! No benchmark output
      OU_BENCHLOG = 0

    end if

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_done ()

!<description>
  ! Closes all log files, cleans up output system.
!</description>

!</subroutine>
  
    ! Close all open channels.
    if (OU_LOG .gt. 0) close(OU_LOG)
    if (OU_ERRORLOG .gt. 0) close(OU_ERRORLOG)
    if (OU_BENCHLOG .gt. 0) close(OU_BENCHLOG)
    
  end subroutine

!************************************************************************************

!<subroutine>

  function output_reformatMsg (smessage, coutputClass, ssubroutine) result (s)

!<description>
  ! Reformats an output message according to an output class.
!</description>

!<input>
  ! The message to reformat.
  character(LEN=*), intent(IN) :: smessage

  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  integer, intent(IN), optional :: coutputClass

  ! OPTIONAL: Name of a subroutine to include in the message; may be "".
  character(LEN=*), intent(IN), optional :: ssubroutine
!</input>

!<result>
  ! Reformatted string
!</result>

  character(LEN=len(smessage)+20+SYS_NAMELEN) :: s

!</subroutine>

    logical :: bsub
    integer :: ioc
    
    bsub = .false.
    if (present (ssubroutine)) then
      bsub = (ssubroutine .ne. '')
    end if
    
    ioc = OU_CLASS_MSG
    if (present(coutputClass)) ioc = coutputClass
      
    if (.not. bsub) then

      select case (ioc)
      case (OU_CLASS_TRACE1)
        s = '*** '//smessage
      case (OU_CLASS_TRACE2)
        s = '***** '//smessage
      case (OU_CLASS_TRACE3)
        s = '******* '//smessage
      case (OU_CLASS_SYSTEM)
        s = '* System: '//smessage
      case (OU_CLASS_ERROR)
        s = '* Error: '//smessage
      case (OU_CLASS_WARNING)
        s = '* Warning: '//smessage
      case DEFAULT
        s = smessage
      end select
    
    else

      select case (ioc)
      case (OU_CLASS_TRACE1)
        s = '*** '//trim(ssubroutine)//': '//smessage
      case (OU_CLASS_TRACE2)
        s = '***** '//trim(ssubroutine)//': '//smessage
      case (OU_CLASS_TRACE3)
        s = '******* '//trim(ssubroutine)//': '//smessage
      case (OU_CLASS_SYSTEM)
        s = '* System ('//trim(ssubroutine)//'): '//smessage
      case (OU_CLASS_ERROR)
        s = '* Error ('//trim(ssubroutine)//'): '//smessage
      case (OU_CLASS_WARNING)
        s = '* Warning ('//trim(ssubroutine)//'): '//smessage
      case DEFAULT
        s = smessage
      end select

    end if

  end function


!************************************************************************************

!<subroutine>

  subroutine output_line_std (smessage, &
                              coutputClass, coutputMode, ssubroutine, &
                              bnolinebreak, bnotrim)

!<description>
  ! Writes an output message to the terminal, log file or error log file,
  ! depending on the input parameters.
  ! smessage is the message to be written out. 
  ! coutputMode decides (if given) about whether the output if written to file,
  ! terminal or both.
  ! coutputClass classifies (if given) the message as standard message, trace
  ! or error message.
!</description>

!<input>
  ! The message to be written out.
  character(LEN=*), intent(IN) :: smessage
  
  ! OPTIONAL: Output mode. One of the OU_MODE_xxxx constants. If not specified,
  ! OU_MODE_STD is assumed.
  integer(I32), intent(IN), optional :: coutputMode

  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  integer, intent(IN), optional :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  character(LEN=*), intent(IN), optional :: ssubroutine
  
  ! OPTIONAL: When specifying bnolinebreak=TRUE, the output routine will
  ! not perform a line break after printing.
  logical, intent(IN), optional :: bnolinebreak

  ! OPTIONAL: When specifying bnotrim=TRUE, the output routine will
  ! not trim the string when printing. This does only work for 
  ! coutputClass=OU_CLASS_MSG and coutputClass=OU_CLASS_ERROR
  ! (or if coutputClass is not specified)!
  logical, intent(IN), optional :: bnotrim
!</input>

!</subroutine>

  ! local variables
  integer(I32) :: coMode
  integer :: coClass, iofChannel, iotChannel
  logical :: bntrim, bnnewline
  character(LEN=len(smessage)+20+SYS_NAMELEN) :: smsg
  
    ! Get the actual parameters.
    
    coMode = OU_MODE_STD
    coClass = OU_CLASS_MSG
    bntrim = .false.
    bnnewline = .false.
    
    if (present(coutputMode))  coMode = coutputMode
    if (present(coutputClass)) coClass = coutputClass
    if (present(bnotrim))      bntrim = bnotrim
    if (present(bnolinebreak)) bnnewline = bnolinebreak
    
    ! Get the file and terminal output channel
    iotChannel = OU_TERMINAL
    if ((coClass .eq. OU_CLASS_ERROR) .or. &
        (coClass .eq. OU_CLASS_WARNING)) iotChannel = OU_ERROR

    iofChannel = OU_LOG
    if ((coClass .eq. OU_CLASS_ERROR) .or. &
        (coClass .eq. OU_CLASS_WARNING)) iofChannel = OU_ERRORLOG
    
    if (bntrim .and. &
        ((coClass .eq. OU_CLASS_MSG) .or. &
         (coClass .eq. OU_CLASS_WARNING) .or. &
         (coClass .eq. OU_CLASS_ERROR)) ) then
      
      ! Where to write the message to?
      if ((iand(coMode,OU_MODE_TERM) .ne. 0) .and. (iotChannel .gt. 0)) then
        if (bnnewline) then
          write (iotChannel,'(A)',ADVANCE='NO') smessage
        else
          write (iotChannel,'(A)') smessage
        end if
      end if
      
      ! Output to log file?
      if ((iand(coMode,OU_MODE_LOG) .ne. 0) .and. (iofChannel .gt. 0)) then
        if (bnnewline) then
          write (iofChannel,'(A)',ADVANCE='NO') smessage
        else
          write (iofChannel,'(A)') smessage
        end if
      end if

      ! Output to the benchmark log file?
      if ((iand(coMode,OU_MODE_BENCHLOG) .ne. 0) .and. (OU_BENCHLOG .gt. 0)) then
        if (bnnewline) then
          write (OU_BENCHLOG,'(A)',ADVANCE='NO') smessage
        else
          write (OU_BENCHLOG,'(A)') smessage
        end if
      end if
      
    else
    
      ! Build the actual error message
      smsg = output_reformatMsg (smessage, coutputClass, ssubroutine) 
      
      ! Where to write the message to?
      if ((iand(coMode,OU_MODE_TERM) .ne. 0) .and. (iotChannel .gt. 0)) then
        if (bnnewline) then
          write (iotChannel,'(A)',ADVANCE='NO') trim(smsg)
        else
          write (iotChannel,'(A)') trim(smsg)
        end if
      end if
      
      ! Output to log file?
      if ((iand(coMode,OU_MODE_LOG) .ne. 0) .and. (iofChannel .gt. 0)) then
        if (bnnewline) then
          write (iofChannel,'(A)',ADVANCE='NO') trim(smsg)
        else
          write (iofChannel,'(A)') trim(smsg)
        end if
      end if

      ! Output to benchmark log file?
      if ((iand(coMode,OU_MODE_BENCHLOG) .ne. 0) .and. (OU_BENCHLOG .gt. 0)) then
        if (bnnewline) then
          write (OU_BENCHLOG,'(A)',ADVANCE='NO') trim(smsg)
        else
          write (OU_BENCHLOG,'(A)') trim(smsg)
        end if
      end if
    end if

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_line_feast (coutputClass, ssubroutine, smsg)

!<description>
  ! Writes an output message to the terminal, log file or error log file,
  ! depending on the input parameters.
  ! smsg is the message to be written out. 
  ! coutputMode decides (if given) about whether the output if written to file,
  ! terminal or both.
  ! coutputClass classifies (if given) the message as standard message, trace
  ! or error message.
  !
  ! Compatibility function for FEAST routines.
!</description>

!<input>
  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  integer, intent(IN) :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  character(LEN=*), intent(IN) :: ssubroutine
  
  ! The message to be written out.
  character(LEN=*), intent(IN) :: smsg
  
!</input>

!</subroutine>

    ! REMARK: Don't rename 'smsg' to 'smessage' -- this does not comply
    ! with the Fortran standard and will give you an ambigious interface
    ! in output_line as Fortran cannot distinguish between output_line_std
    ! and output_line_feast anymore in that case! (as all parameter names
    ! are the same).

    call output_line_std (smsg, coutputClass, OU_MODE_STD, ssubroutine)

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_lbrk (coutputClass, coutputMode, ssubroutine, nlbrk)

!<description>
  ! Writes a line break to the terminal, log file or error log file,
  ! depending on the input parameters.
  ! coutputMode decides (if given) about whether the output if written to file,
  ! terminal or both.
  ! coutputClass classifies (if given) the message as standard message, trace
  ! or error message.
!</description>

!<input>
  ! OPTIONAL: Output mode. One of the OU_MODE_xxxx constants. If not specified,
  ! OU_MODE_STD is assumed.
  integer(I32), intent(IN), optional :: coutputMode

  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  integer, intent(IN), optional :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  character(LEN=*), intent(IN), optional :: ssubroutine
  
  ! OPTIONAL: Number of linebreaks
  integer, intent(IN), optional :: nlbrk
!</input>

!</subroutine>

    ! local variables
    integer :: ilbrk

    if (present(nlbrk)) then

      do ilbrk = 1, nlbrk
        call output_line_std ('', coutputClass, coutputMode, ssubroutine)
      end do

    else

      call output_line_std ('', coutputClass, coutputMode, ssubroutine)

    end if

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_separator (csepType, coutputClass, coutputMode, ssubroutine)

!<description>
  ! Writes a separator line to the terminal, log file or error log file,
  ! depending on the input parameters. The parameter csepType decides
  ! on the type of separator line.
  ! coutputMode decides (if given) about whether the output if written to file,
  ! terminal or both.
  ! coutputClass classifies (if given) the message as standard message, trace
  ! or error message.
!</description>

!<input>
  ! Type of the separator line. One of the OU_SEP_xxxx constants
  ! (OU_SEP_MINUS, OU_SEP_STAR, OU_SEP_EQUAL,...).
  integer, intent(IN) :: csepType

  ! OPTIONAL: Output mode. One of the OU_MODE_xxxx constants. If not specified,
  ! OU_MODE_STD is assumed.
  integer(I32), intent(IN), optional :: coutputMode

  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  integer, intent(IN), optional :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  character(LEN=*), intent(IN), optional :: ssubroutine
  
  ! OPTIONAL: Name
!</input>

!</subroutine>

    ! local variables
    character(LEN=max(1,OU_LINE_LENGTH-1)) :: cstr
    character(LEN=20) :: saux

    ! Create the separator line
    select case (csepType)
    case (OU_SEP_MINUS)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''-''))'
    case (OU_SEP_PLUS)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''+''))'
    case (OU_SEP_STAR)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''*''))'
    case (OU_SEP_EQUAL)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''=''))'
    case (OU_SEP_DOLLAR)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''$''))'
    case (OU_SEP_AT)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''@''))'
    case (OU_SEP_TILDE)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''~''))'
    case (OU_SEP_AMPAND)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''&''))'
    case (OU_SEP_PERC)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''%''))'
    case (OU_SEP_HASH)
      write (saux,'(A,I3,A)') '(',len(cstr),'(''#''))'
    case DEFAULT
      write (unit=*, fmt=*) 'output_separator: Unknown separator type: ',&
                            csepType
      call sys_halt()
    end select
    
    write (cstr,saux)

    call output_line_std (cstr, coutputClass, coutputMode, ssubroutine)

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_simple (ioutputLevel, smsg, coutputClass, ssubroutine)

!<description>
  ! Writes an output message to the terminal, log file or error log file,
  ! depending on the input parameters.
  ! smsg is the message to be written out. 
  ! ioutputLevel decides on where the message is written to, to the terminal
  ! and/or the log file.
  ! coutputClass classifies (if given) the message as standard message, trace
  ! or error message.
!</description>

!<input>
  ! Output level. DEcides on where to write the message to.
  ! <=0: No output,
  !  =1: write to log file.
  ! >=2: write to log file and terminal.
  integer, intent(IN) :: ioutputLevel

  ! OPTIONAL: The message to be written out.
  character(LEN=*), intent(IN), optional :: smsg
  
  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  integer, intent(IN), optional :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  character(LEN=*), intent(IN), optional :: ssubroutine
  
!</input>

!</subroutine>

  integer(I32) :: coMode
  
    select case (ioutputLevel)
    case (:0)
      coMode = 0
    case (1)
      coMode = OU_MODE_LOG
    case (2:)
      coMode = OU_MODE_LOG+OU_MODE_TERM
    end select

    if (present(smsg)) then
      call output_line_std (smsg, coutputClass, coMode, ssubroutine)
    else
      call output_line_std ('', coutputClass, coMode, ssubroutine)
    end if

  end subroutine

!************************************************************************************

!<subroutine>

  subroutine output_simple_sep (ioutputLevel, csepType, coutputClass, ssubroutine)

!<description>
  ! Writes a separator line to the terminal, log file or error log file,
  ! depending on the input parameters. The parameter csepType decides
  ! on the type of separator line.
  ! ioutputLevel decides on where the message is written to, to the terminal
  ! and/or the log file.
  ! coutputClass classifies (if given) the message as standard message, trace
  ! or error message.
!</description>

!<input>
  ! Output level. DEcides on where to write the message to.
  ! <=0: No output,
  !  =1: write to log file.
  ! >=2: write to log file and terminal.
  integer, intent(IN) :: ioutputLevel

  ! Type of the separator line. One of the OU_SEP_xxxx constants
  ! (OU_SEP_MINUS, OU_SEP_STAR or OU_SEP_EQUAL).
  integer, intent(IN) :: csepType
  
  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  integer, intent(IN), optional :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  character(LEN=*), intent(IN), optional :: ssubroutine
  
!</input>

!</subroutine>

  integer(I32) :: coMode
  
    select case (ioutputLevel)
    case (:0)
      coMode = 0
    case (1)
      coMode = OU_MODE_LOG
    case (2:)
      coMode = OU_MODE_LOG+OU_MODE_TERM
    end select

    call output_separator (csepType, coutputClass, coMode, ssubroutine)

  end subroutine
  
end module
