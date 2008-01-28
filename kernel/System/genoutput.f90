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
!#     f) CALL output_line_sep (OU_SEP_MINUS)
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
!# </purpose>
!##############################################################################

MODULE genoutput

  USE fsystem

  IMPLICIT NONE

!<constants>

!<constantblock description="Output mode. Decides on where the output is written to.">

  ! Output mode: Write to main log file / error log file (depending on whether
  ! a message is considered as an error or not by OU_CLASS_ERROR/OU_CLASS_WARNING)
  INTEGER, PARAMETER :: OU_MODE_LOG      = 2**0
  
  ! Output mode: Write to terminal
  INTEGER, PARAMETER :: OU_MODE_TERM     = 2**1

  ! Output mode: Write to both, log file and terminal
  INTEGER, PARAMETER :: OU_MODE_STD      = OU_MODE_LOG+OU_MODE_TERM
  
!</constantblock>
  
!<constantblock description="Output classification. Prints an additional classification string\\
!                            when writing a string to a log file / terminal">

  ! Output classification: Standard message
  INTEGER, PARAMETER :: OU_CLASS_MSG      = 0
  
  ! Output classification: Trace information, level 1
  INTEGER, PARAMETER :: OU_CLASS_TRACE1   = 1

  ! Output classification: Trace information, level 2
  INTEGER, PARAMETER :: OU_CLASS_TRACE2   = 2

  ! Output classification: Trace information, level 3
  INTEGER, PARAMETER :: OU_CLASS_TRACE3   = 3

  ! Output classification: System/Timer message
  INTEGER, PARAMETER :: OU_CLASS_SYSTEM   = 4
  
  ! Output classification: Error message.
  INTEGER, PARAMETER :: OU_CLASS_ERROR    = 5

  ! Output classification: Warning message.
  INTEGER, PARAMETER :: OU_CLASS_WARNING  = 6

!</constantblock>

!<constantblock description="Type of separator line">

  ! Separator line: MINUS character
  INTEGER, PARAMETER :: OU_SEP_MINUS  = 0

  ! Separator line: STAR character 
  INTEGER, PARAMETER :: OU_SEP_STAR   = 1

  ! Separator line: EQUAL character
  INTEGER, PARAMETER :: OU_SEP_EQUAL  = 2

  ! Separator line: DOLLAR character
  INTEGER, PARAMETER :: OU_SEP_DOLLAR = 3

  ! Separator line: @ character
  INTEGER, PARAMETER :: OU_SEP_AT     = 4

  ! Separator line: + character
  INTEGER, PARAMETER :: OU_SEP_PLUS   = 5

!</constantblock>

!<constantblock>

  ! Length of a line on the terminal / in the log file.
  ! Standard value = 80 characters
  INTEGER :: OU_LINE_LENGTH         = 80

  ! Global device number for terminal output
  INTEGER :: OU_TERMINAL            = 6

  ! Global device number for terminal input
  INTEGER :: IN_TERMINAL            = 5

  ! Global device number for errors on terminal
  INTEGER :: OU_ERROR               = 6
  
  ! Global device number for log file. The log file is opened in output_init.
  ! <=0: no log file output.
  INTEGER :: OU_LOG                 = 0
  
  ! Global device number for error log file. The error log file is opened 
  ! in output_init and usually coincides with OU_LOG.
  ! <=0: write error messages to standard log file.
  INTEGER :: OU_ERRORLOG            = 0

!</constantblock>

!</constants>  

  INTERFACE output_line
    MODULE PROCEDURE output_line_std
    MODULE PROCEDURE output_line_feast
  END INTERFACE
  
  INTERFACE output_init
    MODULE PROCEDURE output_init_simple
    MODULE PROCEDURE output_init_logfile
    MODULE PROCEDURE output_init_standard
  END INTERFACE

CONTAINS

!************************************************************************

!<subroutine>
  SUBROUTINE output_openLogfile(sfilename, iunit)

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
    CHARACTER(*), INTENT(IN) :: sfilename

  !</input>

  !<output>
    !unit of the opened file
    INTEGER, INTENT(OUT) :: iunit
  !</output>
  
!</subroutine>

    INTEGER :: istatus ! status variable for opening procedure

    IF (trim(sfilename) .eq. "") THEN 
      WRITE (*,'(A)') 'Error: io_openFileForWriting. sfilename undefined!'
      RETURN
    ENDIF

    iunit = sys_getFreeUnit()
    IF (iunit .eq. -1) THEN
      WRITE (*,'(A)') 'Error: io_openFileForWriting. No output channel available!'
      RETURN
    ENDIF

    OPEN(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
          action="write")

    IF (istatus .ne. 0) THEN
      WRITE(unit=*,fmt=*) 'Error: io_openFileForWriting. &
                          &Error while opening file "', trim(sfilename), '". ***'
      iunit = -1
    ENDIF

  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_init_simple ()

!<description>
  ! Initialises the output system. The output system is configured to show all
  ! messages on the terminal, no log file is opened.
!</description>
  
!</subroutine>

    CALL output_init_standard ("","")

  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_init_logfile (slogFilename)

!<description>
  ! Initialises the output system. 
  ! If sfilename is given, a log file with that name is opened and all log and error
  ! messages are written to it. If not given, all output is directed to the 
  ! terminal only.
!</description>
  
!<input>

  ! Name of a log file for standard messages. If "" is specified, 
  ! output messages are written to the standard output.
  CHARACTER(LEN=*), INTENT(IN) :: slogFilename
  
!</input>

!</subroutine>

    CALL output_init_standard (slogFilename,"")

  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_init_standard (slogFilename,serrorFilename)

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
  CHARACTER(LEN=*), INTENT(IN) :: slogFilename
  
  ! Name of an log file for error messages. If ""
  ! is specified, errors are written to the standard log file. The name of 
  ! the file may also coincide with slogFilename.
  CHARACTER(LEN=*), INTENT(IN) :: serrorFilename
  
!</input>

!</subroutine>

  ! Both filenames given?
  IF ((slogFilename .NE. "") .AND. (serrorFilename .NE. "")) THEN
  
    ! Both filenames the same=
    IF (slogFilename .EQ. serrorFilename) THEN
    
      ! Only open one log file for both.
      CALL output_openLogfile(slogFilename, OU_LOG)
      OU_ERRORLOG = OU_LOG
    
    ELSE
    
      ! Open two separate log files
      CALL output_openLogfile(slogFilename, OU_LOG)
      CALL output_openLogfile(serrorFilename, OU_ERRORLOG)
    
    END IF
  
  ELSE 
  
    ! Close previously opened log files.
    CALL output_done ()    
    
    IF (slogFilename .NE. "") THEN
    
      ! Open a log file
      CALL output_openLogfile(slogFilename, OU_LOG)
    
    END IF
  
    IF (serrorFilename .NE. "") THEN
    
      ! Open an error log file
      CALL output_openLogfile(serrorFilename, OU_ERRORLOG)
      
    ELSE
    
      ! Write error messages to standard log file.
      OU_ERRORLOG = OU_LOG
    
    END IF
  
  END IF
  
  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_done ()

!<description>
  ! Closes all log files, cleans up output system.
!</description>

!</subroutine>
  
    ! Close all open channels.
    IF (OU_LOG .GT. 0) CLOSE(OU_LOG)
    IF (OU_ERRORLOG .GT. 0) CLOSE(OU_ERRORLOG)
    
  END SUBROUTINE

!************************************************************************************

!<subroutine>

  FUNCTION output_reformatMsg (smessage, coutputClass, ssubroutine) RESULT (s)

!<description>
  ! Reformats an output message according to an output class.
!</description>

!<input>
  ! The message to reformat.
  CHARACTER(LEN=*), INTENT(IN) :: smessage

  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: coutputClass

  ! OPTIONAL: Name of a subroutine to include in the message; may be "".
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssubroutine
!</input>

!<result>
  ! Reformatted string
!</result>

  CHARACTER(LEN=LEN(smessage)+20+SYS_NAMELEN) :: s

!</subroutine>

    LOGICAL :: bsub
    INTEGER :: ioc
    
    bsub = .FALSE.
    IF (PRESENT (ssubroutine)) THEN
      bsub = (ssubroutine .NE. '')
    END IF
    
    ioc = OU_CLASS_MSG
    IF (PRESENT(coutputClass)) ioc = coutputClass
      
    IF (.NOT. bsub) THEN

      SELECT CASE (ioc)
      CASE (OU_CLASS_TRACE1)
        s = '*** '//smessage
      CASE (OU_CLASS_TRACE2)
        s = '***** '//smessage
      CASE (OU_CLASS_TRACE3)
        s = '******* '//smessage
      CASE (OU_CLASS_SYSTEM)
        s = '* System: '//smessage
      CASE (OU_CLASS_ERROR)
        s = '* Error: '//smessage
      CASE (OU_CLASS_WARNING)
        s = '* Warning: '//smessage
      CASE DEFAULT
        s = smessage
      END SELECT
    
    ELSE

      SELECT CASE (ioc)
      CASE (OU_CLASS_TRACE1)
        s = '*** '//trim(ssubroutine)//': '//smessage
      CASE (OU_CLASS_TRACE2)
        s = '***** '//trim(ssubroutine)//': '//smessage
      CASE (OU_CLASS_TRACE3)
        s = '******* '//trim(ssubroutine)//': '//smessage
      CASE (OU_CLASS_SYSTEM)
        s = '* System ('//trim(ssubroutine)//'): '//smessage
      CASE (OU_CLASS_ERROR)
        s = '* Error ('//trim(ssubroutine)//'): '//smessage
      CASE (OU_CLASS_WARNING)
        s = '* Warning ('//trim(ssubroutine)//'): '//smessage
      CASE DEFAULT
        s = smessage
      END SELECT

    END IF

  END FUNCTION


!************************************************************************************

!<subroutine>

  SUBROUTINE output_line_std (smessage, &
                              coutputClass, coutputMode, ssubroutine, &
                              bnolinebreak,bnotrim)

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
  CHARACTER(LEN=*), INTENT(IN) :: smessage
  
  ! OPTIONAL: Output mode. One of the OU_MODE_xxxx constants. If not specified,
  ! OU_MODE_STD is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: coutputMode

  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssubroutine
  
  ! OPTIONAL: When specifying bnolinebreak=TRUE, the output routine will
  ! not perform a line break after printing.
  LOGICAL, INTENT(IN), OPTIONAL :: bnolinebreak

  ! OPTIONAL: When specifying bnotrim=TRUE, the output routine will
  ! not trim the string when printing. This does only work for 
  ! coutputClass=OU_CLASS_MSG and coutputClass=OU_CLASS_ERROR
  ! (or if coutputClass is not specified)!
  LOGICAL, INTENT(IN), OPTIONAL :: bnotrim
!</input>

!</subroutine>

  ! local variables
  INTEGER :: coMode, coClass, iofChannel, iotChannel
  LOGICAL :: bntrim, bnnewline
  CHARACTER(LEN=LEN(smessage)+20+SYS_NAMELEN) :: smsg
  
    ! Get the actual parameters.
    
    coMode = OU_MODE_STD
    coClass = OU_CLASS_MSG
    bntrim = .FALSE.
    bnnewline = .FALSE.
    
    IF (PRESENT(coutputMode))  coMode = coutputMode
    IF (PRESENT(coutputClass)) coClass = coutputClass
    IF (PRESENT(bnotrim))      bntrim = bnotrim
    IF (PRESENT(bnolinebreak)) bnnewline = bnolinebreak
    
    ! Get the file and terminal output channel
    iotChannel = OU_TERMINAL
    IF ((coClass .EQ. OU_CLASS_ERROR) .OR. &
        (coClass .EQ. OU_CLASS_WARNING)) iofChannel = OU_ERROR

    iofChannel = OU_LOG
    IF ((coClass .EQ. OU_CLASS_ERROR) .OR. &
        (coClass .EQ. OU_CLASS_WARNING)) iofChannel = OU_ERRORLOG
    
    IF (bntrim .AND. &
        ((coClass .EQ. OU_CLASS_MSG) .OR. &
         (coClass .EQ. OU_CLASS_WARNING) .OR. &
         (coClass .EQ. OU_CLASS_ERROR)) ) THEN
         
      ! Where to write the message to?
      IF ((IAND(coMode,OU_MODE_TERM) .NE. 0) .AND. (iotChannel .GT. 0)) THEN
        IF (bnnewline) THEN
          WRITE (iotChannel,'(A)',ADVANCE='NO') smessage
        ELSE
          WRITE (iotChannel,'(A)') smessage
        END IF
      END IF
      
      ! Output to log file?
      IF ((IAND(coMode,OU_MODE_LOG) .NE. 0) .AND. (iofChannel .GT. 0)) THEN
        IF (bnnewline) THEN
          WRITE (iofChannel,'(A)',ADVANCE='NO') smessage
        ELSE
          WRITE (iofChannel,'(A)') smessage
        END IF
      END IF
      
    ELSE
    
      ! Build the actual error message
      smsg = output_reformatMsg (smessage, coutputClass, ssubroutine) 
      
      ! Where to write the message to?
      IF ((IAND(coMode,OU_MODE_TERM) .NE. 0) .AND. (iotChannel .GT. 0)) THEN
        IF (bnnewline) THEN
          WRITE (iotChannel,'(A)',ADVANCE='NO') TRIM(smsg)
        ELSE
          WRITE (iotChannel,'(A)') TRIM(smsg)
        END IF
      END IF
      
      ! Output to log file?
      IF ((IAND(coMode,OU_MODE_LOG) .NE. 0) .AND. (iofChannel .GT. 0)) THEN
        IF (bnnewline) THEN
          WRITE (iofChannel,'(A)',ADVANCE='NO') TRIM(smsg)
        ELSE
          WRITE (iofChannel,'(A)') TRIM(smsg)
        END IF
      END IF
    END IF

  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_line_feast (coutputClass, ssubroutine, smsg)

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
  INTEGER, INTENT(IN) :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  CHARACTER(LEN=*), INTENT(IN) :: ssubroutine
  
  ! The message to be written out.
  CHARACTER(LEN=*), INTENT(IN) :: smsg
  
!</input>

!</subroutine>

    ! REMARK: Don't rename 'smsg' to 'smessage' -- this does not comply
    ! with the Fortran standard and will give you an ambigious interface
    ! in output_line as Fortran cannot distinguish between output_line_std
    ! and output_line_feast anymore in that case! (as all parameter names
    ! are the same).

    CALL output_line_std (smsg, coutputClass, OU_MODE_STD, ssubroutine)

  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_lbrk (coutputClass, coutputMode, ssubroutine)

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
  INTEGER, INTENT(IN), OPTIONAL :: coutputMode

  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssubroutine
  
  ! OPTIONAL: Name
!</input>

!</subroutine>

    CALL output_line_std ('', coutputClass, coutputMode, ssubroutine)

  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_separator (csepType, coutputClass, coutputMode, ssubroutine)

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
  INTEGER, INTENT(IN) :: csepType

  ! OPTIONAL: Output mode. One of the OU_MODE_xxxx constants. If not specified,
  ! OU_MODE_STD is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: coutputMode

  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssubroutine
  
  ! OPTIONAL: Name
!</input>

!</subroutine>

    ! local variables
    CHARACTER(LEN=MAX(1,OU_LINE_LENGTH-1)) :: cstr
    CHARACTER(LEN=20) :: saux

    ! Create the separator line
    SELECT CASE (csepType)
    CASE (OU_SEP_MINUS)
      WRITE (saux,'(A,I3,A)') '(',LEN(cstr),'(''-''))'
    CASE (OU_SEP_PLUS)
      WRITE (saux,'(A,I3,A)') '(',LEN(cstr),'(''+''))'
    CASE (OU_SEP_STAR)
      WRITE (saux,'(A,I3,A)') '(',LEN(cstr),'(''*''))'
    CASE (OU_SEP_EQUAL)
      WRITE (saux,'(A,I3,A)') '(',LEN(cstr),'(''=''))'
    CASE (OU_SEP_DOLLAR)
      WRITE (saux,'(A,I3,A)') '(',LEN(cstr),'(''$''))'
    CASE (OU_SEP_AT)
      WRITE (saux,'(A,I3,A)') '(',LEN(cstr),'(''@''))'
    CASE DEFAULT
      PRINT *,'output_separator: Unknown separator type: ',csepType
      CALL sys_halt()
    END SELECT
    
    WRITE (cstr,saux)

    CALL output_line_std (cstr, coutputClass, coutputMode, ssubroutine)

  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_simple (ioutputLevel, smsg, coutputClass, ssubroutine)

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
  INTEGER, INTENT(IN) :: ioutputLevel

  ! OPTIONAL: The message to be written out.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: smsg
  
  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssubroutine
  
!</input>

!</subroutine>

  INTEGER :: coMode
  
    SELECT CASE (ioutputLevel)
    CASE (:0)
      coMode = 0
    CASE (1)
      coMode = OU_MODE_LOG
    CASE (2:)
      coMode = OU_MODE_LOG+OU_MODE_TERM
    END SELECT

    IF (PRESENT(smsg)) THEN
      CALL output_line_std (smsg, coMode, coutputClass, ssubroutine)
    ELSE
      CALL output_line_std ('', coMode, coutputClass, ssubroutine)
    END IF

  END SUBROUTINE

!************************************************************************************

!<subroutine>

  SUBROUTINE output_simple_sep (ioutputLevel, csepType, coutputClass, ssubroutine)

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
  INTEGER, INTENT(IN) :: ioutputLevel

  ! Type of the separator line. One of the OU_SEP_xxxx constants
  ! (OU_SEP_MINUS, OU_SEP_STAR or OU_SEP_EQUAL).
  INTEGER, INTENT(IN) :: csepType
  
  ! OPTIONAL: Output classification. One of the OU_CLASS_xxxx constants.
  ! If not specified, OU_CLASS_MSG is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: coutputClass
  
  ! OPTIONAL: Name of the subroutine that calls this function
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssubroutine
  
!</input>

!</subroutine>

  INTEGER :: coMode
  
    SELECT CASE (ioutputLevel)
    CASE (:0)
      coMode = 0
    CASE (1)
      coMode = OU_MODE_LOG
    CASE (2:)
      coMode = OU_MODE_LOG+OU_MODE_TERM
    END SELECT

    CALL output_separator (csepType, coMode, coutputClass, ssubroutine)

  END SUBROUTINE
  
END MODULE
