!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2scriptfile </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a simple parser for text commands from a script file
!# or from the terminal. Once invoked, the scr_readScript and/or 
!# scr_readTerminal routines will parse commands entered by a script file
!# or from the terminal. This allows on-line modification during the
!# execution of the program.
!# 
!# The following routines can be found here:
!#
!# 1.) scr_readScript
!#     -> Reads in a text file and interprets every line as a command.
!#        Returns if all commands are processed.
!#
!# 2.) scr_readTerminal
!#     -> Starts text input from the terminal; returns if the user
!#        enters a CONTINUE command on the terminal.
!#
!# 3.) scr_parseLine
!#     -> This is the central parser subroutine. This routine interprets
!#        a line that comes from the terminal or from a text file and
!#        processes this command.
!#
!# </purpose>
!##############################################################################

MODULE cc2dminim2scriptfile

  USE fsystem
  USE storage
  
  USE cc2dmediumm2basic
  
  IMPLICIT NONE

CONTAINS

! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE scr_readScript (sfilename,ccleardelete,rproblem)

!<description>
  ! Reads a text file sfilename line by line and calls scr_parseLine
  ! for each of the lines. 
!</description>

!<input>

  ! Name of the file to be read
  CHARACTER(LEN=*), INTENT(IN) :: sfilename

  ! Whether to delete or clear the file after reading
  ! =0: don't do anything.
  ! =1: clear the file.
  ! =2: delete the file.
  INTEGER, INTENT(IN) :: ccleardelete

!</input>

!<inputoutput>

  ! Problem structure; passed to scr_parseLine.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

!</inputoutput>

!</subroutine>

    ! Open the file -- exclusively for us!
    ! This is somehow like io_openFileForReading but cancels if
    ! the file cannot be exclusively opened.

    LOGICAL :: bexists, bcontinue
    INTEGER :: istatus,iunit,ilinelen,istat
    CHARACTER(LEN=4096) :: sdata

    IF (trim(sfilename) .eq. '') THEN
      CALL output_line ('No file name!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'scr_readScript')
      CALL sys_halt()
    ENDIF

    iunit = sys_getFreeUnit()
    IF (iunit .eq. -1) THEN
      CALL output_line ('No free unit number!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'scr_readScript')
      RETURN
    ENDIF

    ! Check if the file exists. Open it for reading and writing to get
    ! exclusive access.
    INQUIRE(file=trim(sfilename), exist=bexists)

    IF (bexists) THEN
    
      OPEN(unit=iunit, file=trim(sfilename), iostat=istatus, &
          action="readwrite",form="formatted")
          
      IF (istatus .NE. 0) THEN
        RETURN
      END IF
      
    ELSE
      RETURN
    END IF
    
    ! Read the file, line by line.
    ! Call scr_parseLine for each line.
    bcontinue = .FALSE.
    istat = 0
    DO WHILE (istat .EQ. 0)
      CALL io_readlinefromfile (iunit, sdata, ilinelen, istat)
      IF (ilinelen .NE. 0) THEN
        CALL scr_parseLine (sdata(1:ilinelen),bcontinue,rproblem)
        
        ! Probably leave the loop and continue with the processing.
        IF (bcontinue) EXIT
      END IF
    END DO
    
    ! Close the file, delete it if necessary.
    CLOSE (iunit)
    IF (ccleardelete .EQ. 2) THEN
      CALL io_deleteFile (sfilename)
    ELSE IF (ccleardelete .EQ. 1) THEN
      CALL io_openFileForWriting(sfilename, iunit, SYS_REPLACE)
      CLOSE(iunit)
    END IF

  END SUBROUTINE

! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE scr_readTerminal (rproblem)

!<description>
  ! Reads from the terminal line by line and calls scr_parseLine
  ! for each of the lines. 
!</description>

!<inputoutput>

  ! Problem structure; passed to scr_parseLine.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

!</inputoutput>

!</subroutine>

    INTEGER :: ilinelen,istat
    CHARACTER(LEN=4096) :: sdata
    LOGICAL :: bcontinue

    ! Read the file, line by line.
    ! Call scr_parseLine for each line.
    bcontinue = .FALSE.
    istat = 0
    DO WHILE (.NOT. bcontinue)
      CALL output_line ('> ', bnolinebreak=.TRUE.,bnotrim=.TRUE.)
      CALL io_readlinefromfile (IN_TERMINAL, sdata, ilinelen, istat)
      IF (ilinelen .NE. 0) THEN
        CALL output_line (sdata(1:ilinelen), OU_CLASS_MSG,OU_MODE_LOG)
        CALL scr_parseLine (sdata(1:ilinelen),bcontinue,rproblem)
      END IF
    END DO

  END SUBROUTINE

! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE scr_parseLine (sline,bcontinue,rproblem)

!<description>
  ! Parses a line of a script file and executes the commands in there.
!</description>

!<input>
  ! Line to be parsed.
  CHARACTER(LEN=*), INTENT(IN) :: sline
!</input>

!<inputoutput>
  ! Continuation flag. If this flag is set to TRUE, the routine advises
  ! the caller to stop processing commands and continue with the
  ! calculation.
  LOGICAL, INTENT(INOUT) :: bcontinue

  ! Problem structure.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    CHARACTER(LEN=SYS_NAMELEN) :: scommand,sparam1,sparam2
    
    ! Get the command
    IF (TRIM(sline) .EQ. '') RETURN
    READ (sline,*) scommand
    CALL sys_toupper(scommand)
    
    ! Interpret the command, execute the corresponding subfunction
    !
    ! ----------------------- STANDARD COMMANDS -------------------------------
    IF (scommand .EQ. '') THEN
      RETURN
    END IF
    
    IF (scommand .EQ. 'PAUSE') THEN
      READ *
      RETURN
    END IF
    
    IF (scommand .EQ. 'CONTINUE') THEN
      bcontinue = .TRUE.
      RETURN
    END IF

    IF (scommand .EQ. 'C') THEN
      bcontinue = .TRUE.
      RETURN
    END IF

    IF (scommand .EQ. 'REM') THEN
      ! Ignore the line
      RETURN
    END IF
    
    IF (scommand(1:1) .EQ. '#') THEN
      ! Ignore the line
      RETURN
    END IF
    
    IF (scommand .EQ. 'EXIT') THEN
      ! CONTINUE and EXIT are the same
      bcontinue = .TRUE.
      RETURN
    END IF
    
    IF (scommand .EQ. 'STOP') THEN
      sys_haltmode = SYS_HALT_STOP
      CALL sys_halt()
    END IF
    
    IF (scommand .EQ. 'TERMINAL') THEN
      CALL scr_readTerminal (rproblem)
      RETURN
    END IF

    IF (scommand .EQ. 'SCRIPT') THEN
      ! Get the name of the script file to execute.
      READ (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      READ (sline,*) sparam1,sparam2
      
      ! Execute that script
      CALL scr_readScript (sparam2,0,rproblem)
      RETURN
    END IF
    
    ! ----------------------- APPLICATION-SPECIFIC COMMANDS -------------------
    
    IF (scommand .EQ. 'WRITE_SOLUTION') THEN
      ! Get the name of the file
      READ (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      READ (sline,*) sparam1,sparam2

      CALL write_vector (sparam2,rproblem%rdataOneshot%p_rx)
      RETURN
    END IF

    IF (scommand .EQ. 'READ_SOLUTION') THEN
      ! Get the name of the file
      READ (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      READ (sline,*) sparam1,sparam2

      CALL read_vector(sparam2,rproblem%rdataOneshot%p_rx)
      RETURN
    END IF

    IF (scommand .EQ. 'WRITE_RHS') THEN
      ! Get the name of the file
      READ (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      READ (sline,*) sparam1,sparam2

      CALL write_vector (sparam2,rproblem%rdataOneshot%p_rb)
      RETURN
    END IF

    IF (scommand .EQ. 'READ_RHS') THEN
      ! Get the name of the file
      READ (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      READ (sline,*) sparam1,sparam2

      CALL read_vector(sparam2,rproblem%rdataOneshot%p_rb)
      RETURN
    END IF

    ! -------------------------------------------------------------------------

    CALL output_line('Unknown command!')
    CALL output_lbrk()

  CONTAINS
  
    ! Write a space-time solution to disc
    SUBROUTINE write_vector(sfilename,rxsuper)
      CHARACTER(LEN=*),INTENT(IN) :: sfilename
      TYPE(t_spaceTimeVector), INTENT(INOUT) :: rxsuper
      TYPE(t_vectorBlock) :: rx
      IF (sfilename .EQ. '') THEN
        CALL output_line('Unknown filename!')
        RETURN
      END IF
      CALL sptivec_convertSupervecToVector (rxsuper, rx)
      CALL vecio_writeBlockVectorHR(rx,'vec',.FALSE.,0,sfilename,'(E20.10)')
      CALL lsysbl_releaseVector (rx)
      CALL output_line('Solution written.')
    END SUBROUTINE

    ! Reads a space-time solution from disc
    SUBROUTINE read_vector(sfilename,rxsuper)
      CHARACTER(LEN=*),INTENT(IN) :: sfilename
      CHARACTER(LEN=SYS_NAMELEN) :: sname
      TYPE(t_spaceTimeVector), INTENT(INOUT) :: rxsuper
      TYPE(t_vectorBlock) :: rx
      IF (sfilename .EQ. '') THEN
        CALL output_line('Unknown filename!')
        RETURN
      END IF
      CALL vecio_readBlockVectorHR(rx,sname,.FALSE.,0,sfilename,.TRUE.)
      CALL sptivec_convertVectorToSupervec (rx,rxsuper)
      CALL lsysbl_releaseVector (rx)
      CALL output_line('Solution read.')
    END SUBROUTINE

  END SUBROUTINE

END MODULE
