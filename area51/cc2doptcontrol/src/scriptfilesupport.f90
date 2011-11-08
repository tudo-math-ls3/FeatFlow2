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

module scriptfilesupport

  use fsystem
  use io
  use storage
  
  use basicstructures
  
  implicit none

contains

! ***************************************************************************

!<subroutine>

  recursive subroutine scr_readScript (sfilename,ccleardelete,rproblem)

!<description>
  ! Reads a text file sfilename line by line and calls scr_parseLine
  ! for each of the lines.
!</description>

!<input>

  ! Name of the file to be read
  character(LEN=*), intent(IN) :: sfilename

  ! Whether to delete or clear the file after reading
  ! =0: don't do anything.
  ! =1: clear the file.
  ! =2: delete the file.
  integer, intent(IN) :: ccleardelete

!</input>

!<inputoutput>

  ! Problem structure; passed to scr_parseLine.
  type(t_problem), intent(INOUT) :: rproblem

!</inputoutput>

!</subroutine>

    ! Open the file -- exclusively for us!
    ! This is somehow like io_openFileForReading but cancels if
    ! the file cannot be exclusively opened.

    logical :: bexists, bcontinue
    integer :: istatus,iunit,ilinelen,istat
    character(LEN=4096) :: sdata

    if (trim(sfilename) .eq. '') then
      call output_line ('No file name!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'scr_readScript')
      call sys_halt()
    endif

    iunit = sys_getFreeUnit()
    if (iunit .eq. -1) then
      call output_line ('No free unit number!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'scr_readScript')
      return
    endif

    ! Check if the file exists. Open it for reading and writing to get
    ! exclusive access.
    inquire(file=trim(sfilename), exist=bexists)

    if (bexists) then
    
      open(unit=iunit, file=trim(sfilename), iostat=istatus, &
          action="readwrite",form="formatted")
          
      if (istatus .ne. 0) then
        return
      end if
      
    else
      return
    end if
    
    ! Read the file, line by line.
    ! Call scr_parseLine for each line.
    bcontinue = .false.
    istat = 0
    do while (istat .eq. 0)
      call io_readlinefromfile (iunit, sdata, ilinelen, istat)
      if (ilinelen .ne. 0) then
        call scr_parseLine (sdata(1:ilinelen),bcontinue,rproblem)
        
        ! Probably leave the loop and continue with the processing.
        if (bcontinue) exit
      end if
    end do
    
    ! Close the file, delete it if necessary.
    close (iunit)
    if (ccleardelete .eq. 2) then
      call io_deleteFile (sfilename)
    else if (ccleardelete .eq. 1) then
      call io_openFileForWriting(sfilename, iunit, SYS_REPLACE)
      close(iunit)
    end if

  end subroutine

! ***************************************************************************

!<subroutine>

  recursive subroutine scr_readTerminal (rproblem)

!<description>
  ! Reads from the terminal line by line and calls scr_parseLine
  ! for each of the lines.
!</description>

!<inputoutput>

  ! Problem structure; passed to scr_parseLine.
  type(t_problem), intent(INOUT) :: rproblem

!</inputoutput>

!</subroutine>

    integer :: ilinelen,istat
    character(LEN=4096) :: sdata
    logical :: bcontinue

    ! Read the file, line by line.
    ! Call scr_parseLine for each line.
    bcontinue = .false.
    istat = 0
    do while (.not. bcontinue)
      call output_line ('> ', bnolinebreak=.true.,bnotrim=.true.)
      call io_readlinefromfile (IN_TERMINAL, sdata, ilinelen, istat)
      if (ilinelen .ne. 0) then
        call output_line (sdata(1:ilinelen), OU_CLASS_MSG,OU_MODE_LOG)
        call scr_parseLine (sdata(1:ilinelen),bcontinue,rproblem)
      end if
    end do

  end subroutine

! ***************************************************************************

!<subroutine>

  recursive subroutine scr_parseLine (sline,bcontinue,rproblem)

!<description>
  ! Parses a line of a script file and executes the commands in there.
!</description>

!<input>
  ! Line to be parsed.
  character(LEN=*), intent(IN) :: sline
!</input>

!<inputoutput>
  ! Continuation flag. If this flag is set to TRUE, the routine advises
  ! the caller to stop processing commands and continue with the
  ! calculation.
  logical, intent(INOUT) :: bcontinue

  ! Problem structure.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    character(LEN=SYS_NAMELEN) :: scommand,sparam1,sparam2
    
    ! Get the command
    if (trim(sline) .eq. '') return
    read (sline,*) scommand
    call sys_toupper(scommand)
    
    ! Interpret the command, execute the corresponding subfunction
    !
    ! ----------------------- STANDARD COMMANDS -------------------------------
    if (scommand .eq. '') then
      return
    end if
    
    if (scommand .eq. 'PAUSE') then
      read *
      return
    end if
    
    if (scommand .eq. 'CONTINUE') then
      bcontinue = .true.
      return
    end if

    if (scommand .eq. 'C') then
      bcontinue = .true.
      return
    end if

    if (scommand .eq. 'REM') then
      ! Ignore the line
      return
    end if
    
    if (scommand(1:1) .eq. '#') then
      ! Ignore the line
      return
    end if
    
    if (scommand .eq. 'EXIT') then
      ! CONTINUE and EXIT are the same
      bcontinue = .true.
      return
    end if
    
    if (scommand .eq. 'STOP') then
      sys_haltmode = SYS_HALT_STOP
      call sys_halt()
    end if
    
    if (scommand .eq. 'TERMINAL') then
      call scr_readTerminal (rproblem)
      return
    end if

    if (scommand .eq. 'SCRIPT') then
      ! Get the name of the script file to execute.
      read (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      read (sline,*) sparam1,sparam2
      
      ! Execute that script
      call scr_readScript (sparam2,0,rproblem)
      return
    end if
    
    ! ----------------------- APPLICATION-SPECIFIC COMMANDS -------------------
    
    if (scommand .eq. 'WRITE_SOLUTION') then
      ! Get the name of the file
      read (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      read (sline,*) sparam1,sparam2

      call write_vector (sparam2,rproblem%rdataOneshot%p_rx)
      return
    end if

    if (scommand .eq. 'READ_SOLUTION') then
      ! Get the name of the file
      read (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      read (sline,*) sparam1,sparam2

      call read_vector(sparam2,rproblem%rdataOneshot%p_rx)
      return
    end if

    if (scommand .eq. 'WRITE_RHS') then
      ! Get the name of the file
      read (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      read (sline,*) sparam1,sparam2

      call write_vector (sparam2,rproblem%rdataOneshot%p_rb)
      return
    end if

    if (scommand .eq. 'READ_RHS') then
      ! Get the name of the file
      read (sline,*) scommand,sparam1
      
      ! Cancel away the quotation marks
      read (sline,*) sparam1,sparam2

      call read_vector(sparam2,rproblem%rdataOneshot%p_rb)
      return
    end if

    ! -------------------------------------------------------------------------

    call output_line('Unknown command!')
    call output_lbrk()

  contains
  
    ! Write a space-time solution to disc
    subroutine write_vector(sfilename,rxsuper)
      character(LEN=*),intent(IN) :: sfilename
      type(t_spaceTimeVector), intent(INOUT) :: rxsuper
      type(t_vectorBlock) :: rx
      if (sfilename .eq. '') then
        call output_line('Unknown filename!')
        return
      end if
      call sptivec_convertSupervecToVector (rxsuper, rx)
      call vecio_writeBlockVectorHR(rx,'vec',.false.,0,sfilename,'(E20.10)')
      call lsysbl_releaseVector (rx)
      call output_line('Solution written.')
    end subroutine

    ! Reads a space-time solution from disc
    subroutine read_vector(sfilename,rxsuper)
      character(LEN=*),intent(IN) :: sfilename
      character(LEN=SYS_NAMELEN) :: sname
      type(t_spaceTimeVector), intent(INOUT) :: rxsuper
      type(t_vectorBlock) :: rx
      if (sfilename .eq. '') then
        call output_line('Unknown filename!')
        return
      end if
      call vecio_readBlockVectorHR(rx,sname,.false.,0,sfilename,.true.)
      call sptivec_convertVectorToSupervec (rx,rxsuper)
      call lsysbl_releaseVector (rx)
      call output_line('Solution read.')
    end subroutine

  end subroutine

end module
