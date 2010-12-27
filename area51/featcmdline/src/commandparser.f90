module commandparser

  use iso_vst
  
  use fsystem
  use genoutput
  use storage
  use collection
  use io
  use basicgeometry
  use boundary
  use triangulation
  
  use meshhierarchy
  use element
  use cubature

  use fespacehierarchybase
  use fespacehierarchy
  use spatialdiscretisation
  use multilevelprojection  
  use linearsystemscalar
  use linearsystemblock
  
  use multilevelprojection
  
  use derivatives
  use scalarpde
  use domainintegration
  use bilinearformevaluation
  use stdoperators
  use feevaluation
  use analyticprojection
  use spdiscprojection
  
  use vectorio
  use ucd
  
  use typedsymbol
  use commandparserbase
  use featcommandhandler

  implicit none
  
  private
  
  ! Block of commands.
  type t_commandBlock
  
    ! Number of lines.
    integer :: icount = 0
  
    ! List of lines.
    type(VARYING_STRING), dimension(:), pointer :: p_Rcommands
  end type
  
  public :: cmdprs_init
  public :: cmdprs_done
  public :: cmdprs_parsestream
  public :: cmdprs_parseterminal
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cmdprs_init (rcmdStatus)

  !<description>
    ! Initialises a command line.
  !</description>
  
  !<output>
    ! Command line to initialise.
    type(t_commandstatus), intent(out) :: rcmdStatus
  !</output>

!</subroutine>
  
    call collct_init (rcmdStatus%rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cmdprs_done (rcmdStatus)

  !<description>
    ! Clean up a command line.
  !</description>

  !<inputoutput>
    ! Command line to clean up.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

!</subroutine>
  
    call collct_done (rcmdStatus%rcollection)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cmdprs_initCmdBlock (rcmdBlock)

  !<description>
    ! Initialises a command block.
  !</description>
  
  !<output>
    ! Command block to initialise.
    type(t_commandBlock), intent(out) :: rcmdBlock
  !</output>

!</subroutine>
  
    rcmdBlock%icount = 0
    allocate(rcmdBlock%p_Rcommands(128))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cmdprs_addCommand (rcmdBlock,scommand)

  !<description>
    ! Puts a command to a command block.
  !</description>

  !<inputoutput>
    ! Command block to modify
    type(t_commandBlock), intent(inout) :: rcmdBlock
  !</inputoutput>
  
  !<input>
    ! Command to put into the command block.
    type(VARYING_STRING), intent(in) :: scommand
  !</input>
  
!</subroutine>

    ! local variables
    type(VARYING_STRING), dimension(:), pointer :: p_Rcommands

    ! Reallocate if there is not enough memory.
    if (rcmdBlock%icount .eq. size(rcmdBlock%p_Rcommands)) then
      allocate (p_Rcommands(rcmdBlock%icount+128))
      p_Rcommands(1:rcmdBlock%icount) = rcmdBlock%p_Rcommands(:)
      deallocate(rcmdBlock%p_Rcommands)
      rcmdBlock%p_Rcommands => p_Rcommands
    end if
    
    ! Store the command
    rcmdBlock%icount = rcmdBlock%icount + 1
    rcmdBlock%p_Rcommands(rcmdBlock%icount) = scommand
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cmdprs_doneCmdBlock (rcmdBlock)

  !<description>
    ! Clean up a command block.
  !</description>

  !<inputoutput>
    ! Command block to clean up.
    type(t_commandBlock), intent(inout) :: rcmdBlock
  !</inputoutput>
  
!</subroutine>

    integer :: i
    do i=1,rcmdBlock%icount
      ! Release memory.
      rcmdBlock%p_Rcommands(i) = ""
    end do

    rcmdBlock%icount = 0
    deallocate(rcmdBlock%p_Rcommands)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cmdprs_releaseargs (p_Sargs)
  
  !<description>
    ! Releases memory of a command line array.
  !</description>
  
  !<inputoutput>
    ! Command line to release.
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
  !</inputoutput>

!</subroutine>

    integer :: i
    
    do i=1,size(p_Sargs)
      p_Sargs(i) = ""
    end do
    
    deallocate (p_Sargs)

  end subroutine


  ! ***************************************************************************

  subroutine cmdprs_getparam (Sargs,iparam,sparam,bupcase,bdequote)
  
  !<description>
    ! Returns parameter iparam from a list of parameters.
  !</description>
  
  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
    
    ! Number of the argument.
    integer, intent(in) :: iparam
    
    ! Whether to upper-case the parameter or nor
    logical, intent(in) :: bupcase

    ! Whether to dequote the parameter or not.
    logical, intent(in) :: bdequote
  !</input>
  
  !<output>
    ! Returns the iparam'sth parameter or "" if that one does
    ! not exist.
    character(len=*), intent(out) :: sparam
  !</output>
  
    ! local variables
  
    if (iparam .gt. size(Sargs)) then
      ! Does not exist.
      sparam = ""
      return
    end if
    
    ! Get the parameter
    sparam = adjustl(char(Sargs(iparam)))

    if (bupcase) then
      call sys_toupper(sparam)
    end if
    
    if (bdequote) then
      call sys_dequote(sparam)
    end if

  end subroutine
  
  ! ***************************************************************************

  subroutine cmdprs_splitline (ssource,p_Sargs)
  
  !<description>
    ! Subduvides a line into an array of substrings.
    ! Quotation marks around substrings are removed.
  !</description>
  
  !<input>
    ! A source string.
    type(VARYING_STRING), intent(in) :: ssource
  !</input>
  
  !<output>
    ! A pointer pointing to an array of strings, each string referring to
    ! one parameter in ssource.
    ! The value NULL is returned if there are no arguments
    ! in the string.
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
  !</output>
  
    character(len=(len(ssource)*2+1)) :: stemp,stemp2

    ! Copy the string.
    stemp = ssource
    
    ! Split the string.
    call cmdprs_complexSplit (stemp,stemp2,2,.true.)

    ! Convert to string list.
    call cmdprs_splitlinedirect (stemp2,p_Sargs)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_splitlinedirect (ssource,p_Sargs)
  
  !<description>
    ! Subduvides a line into an array of substrings.
    ! 'Direct' Version, ssource must have been split with cmdprs_complexSplit.
    ! Quotation marks around substrings are removed.
  !</description>
  
  !<input>
    ! A source string.
    character(len=*), intent(in) :: ssource
  !</input>
  
  !<output>
    ! A pointer pointing to an array of strings, each string referring to
    ! one parameter in ssource.
    ! The value NULL is returned if there are no arguments
    ! in the string.
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
  !</output>
  
    integer :: i, icount, istart,iend, imax
    
    nullify(p_Sargs)
    
    ! Stop if there is nothing.
    if (ssource .eq. "") then
      return
    end if
    
    ! Count the number of \0 characters
    icount = 0
    do i = 1,len_trim(ssource)
      if (ssource(i:i) .eq. char(0)) then
        icount = icount + 1
      end if
    end do

    ! Allocate memory.
    allocate (p_Sargs(icount+1))
    
    ! Get the tokens.
    istart = 1
    iend = 0
    i = 0
    imax = len_trim(ssource)
    do while (iend .le. imax)
      ! End of the word
      iend = index (ssource(istart:imax),char(0))
      if (iend .eq. 0) then
        iend = imax+1
      else
        iend = iend+istart-1    ! Correct the position
      end if
      
      ! Save the token
      i = i + 1
      p_Sargs(i) = ssource(istart:iend-1)
      
      ! Next token
      istart = iend + 1
    end do
    
  end subroutine

  ! ***************************************************************************

  recursive subroutine cmdprs_getNextBlock (rcmdblock,iline,istartline,iendline)
  
  !<description>
    ! Parses the lines in a command block and finds a block.
    ! A block is either one single command or a block startng with "{" and ending
    ! with "}".
  !</description>
  
  !<input>
    ! Command block with commands.
    type(t_commandBlock), intent(in) :: rcmdblock
    
    ! Start line from where to parse commands
    integer, intent(in) :: iline
  !</input>
    
  !<output>
    ! Start line in the command block. 
    ! A value =0 indicates an error.
    ! If istartline=0 and iendline>0, the block end was not found.
    integer, intent(out) :: istartline
    
    ! End line in the command block. 
    ! If this equals to istartline, this is a 1-line command block.
    ! Otherwise, the command block starts with "{" and ends with "}".
    integer, intent(out) :: iendline
  !</output>

    integer :: icurrentline,iblockcount
    character(len=32) :: sline
    
    iblockcount = 0
    istartline = 0
    iendline = 0
    
    ! Parse all lines.
    do icurrentline = iline,rcmdblock%icount
    
      ! Get the line
      sline = rcmdblock%p_Rcommands(icurrentline)
      sline = trim(adjustl(sline))
      
      ! Comment? empty line?
      if ((sline .ne. "") .and. (sline(1:1) .ne. "#")) then
      
        ! Command in a block?
        if (iblockcount .gt. 0) then
        
          if (istartline .gt. 0) then
            ! Remember line index for them last line.
            iendline = icurrentline
          end if
          
        end if
        
        ! Block start?
        if (sline .eq. "{") then
        
          ! One more block.
          iblockcount = iblockcount + 1
          
          if (iblockcount .eq. 1) then
            ! First block.
            !
            ! Curretly that's a zero block.
            istartline = icurrentline
            iendline = istartline
          
            ! Skip this line.
            cycle
          end if
          
        else if (sline .eq. "}") then  ! Block end?
        
          ! Decrease block counter
          if (iblockcount .gt. 0) then
            iblockcount = iblockcount - 1
            
            ! Block end?
            if (iblockcount .eq. 0) then
              return
            end if
          else
            ! Error.
            istartline = 0
          end if
          
        end if
      
        ! Command in a block?
        if (iblockcount .eq. 0) then
        
          ! This is a single line command.
          istartline = icurrentline
          iendline = icurrentline
          return
        
        end if        
      
      end if

    end do

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_getNextLine (rcmdblock,iline,istartline,iendline,sline)
  
  !<description>
    ! Increases the line counter to the next line
  !</description>
  
  !<input>
    ! Command block with commands.
    type(t_commandBlock), intent(in) :: rcmdblock

    ! Start line in the command block. 
    integer, intent(in) :: istartline
    
    ! End line in the command block. 
    integer, intent(in) :: iendline
  !</input>
    
  !<inputoutput>
    ! On input: Current line.
    ! On output: Line number of the next valid line or =0 if there is no
    ! more valid line in the block.
    integer, intent(inout) :: iline
  !</inputoutput>
  
  !<output>
    ! The next line. Completely split and dequoted.
    character(len=*) :: sline
  !</output>
  
  !</subroutine>
  
    character(len=len(sline)) :: stemp
    integer :: iprevline
    
    iprevline = iline
    
    do iline = max(istartline,iprevline+1),iendline
    
      ! Get the line and split + dequote it.
      stemp = rcmdblock%p_Rcommands(iline)
      call cmdprs_complexSplit (stemp,sline,1,.false.)
      
      ! Return the first nonempty line
      if (sline .ne. "") return
      
    end do
    
    ! We processed all lines, none found.
    sline = ""
    iline = 0

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  recursive subroutine cmdprs_docommand (rcmdStatus,rcmdblock,inestlevel,&
      scommand,iline,iblockstart,iblockend,bworkflowAllowed,rvalue)
  
  !<description>
    ! Executes a command.
  !</description>
  
  !<inputoutput>
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Input: Number of the line corresponding to scommand.
    ! Output: Last executed line. Program continues after that.
    integer, intent(inout) :: iline
  !</inputoutput>

  !<input>
    ! Command block with commands.
    ! The block may start with "{"; in this case it must end with "}".
    ! the block must at least contain one line.
    type(t_commandBlock), intent(in) :: rcmdblock
    
    ! Command string. Splitted.
    character(len=*), intent(in) :: scommand
    
    ! Start of the commands depending on the command in rcmdBlock.
    integer, intent(in) :: iblockstart

    ! End of the commands depending on the command in rcmdBlock.
    integer, intent(in) :: iblockend
    
    ! Level of nesting
    integer, intent(in) :: inestlevel
    
    ! If set to TRUE, workfloa commands (like FOR, IF, ...) are allowed.
    logical, intent(in) :: bworkflowAllowed
    
    ! Structure encapsuling return values.
    ! Initialised with default values by "intent(out)".
    type(t_symbolValue), intent(out) :: rvalue
  !</input>
  
  !</subroutine>

    integer :: istart, iend, istart2, iend2
    character(len=len_trim(scommand)) :: scommandtrimmed
    logical :: bmatch,btrailexists
    integer :: ierror
    character(len=SYS_NAMELEN) :: ssectionname
    character(len=*), dimension(5), parameter :: spatternBEGIN = &
      (/ "{     ", "      ", "      ", "      ", "      " /)
    character(len=*), dimension(5), parameter :: spatternEND = &
      (/ "}     ", "      ", "      ", "      ", "      " /)
    character(len=*), dimension(5), parameter :: spatternFOR = &
      (/ "FOR   ", "(     ", "*     ", ")     ", "      " /)
    character(len=*), dimension(5), parameter :: spatternWHILE = &
      (/ "WHILE ", "(     ", "*     ", ")     ", "      " /)
    character(len=*), dimension(5), parameter :: spatternDO = &
      (/ "DO    ", "(     ", "*     ", ")     ", "      " /)
    character(len=*), dimension(5), parameter :: spatternIF = &
      (/ "IF    ", "(     ", "*     ", ")     ", "      " /)
    character(len=*), dimension(5), parameter :: spatternINT = &
      (/ "INT   ", "?     ", "      ", "      ", "      " /)
    character(len=*), dimension(5), parameter :: spatternDOUBLE = &
      (/ "DOUBLE", "?     ", "      ", "      ", "      " /)
    character(len=*), dimension(5), parameter :: spatternSTRING = &
      (/ "STRING", "?     ", "      ", "      ", "      " /)
    character(len=*), dimension(5), parameter :: spatternINT2 = &
      (/ "INT   ", "?     ", "=     ", "*     ", "      " /)
    character(len=*), dimension(5), parameter :: spatternDOUBLE2 = &
      (/ "DOUBLE", "?     ", "=     ", "*     ", "      " /)
    character(len=*), dimension(5), parameter :: spatternSTRING2 = &
      (/ "STRING", "?     ", "=     ", "*     ", "      " /)
    
    bmatch = .false.
    
    scommandtrimmed = scommand
    call cmdprs_removeTrail (scommandtrimmed,btrailexists)
    
    ! Name of the section corresponding to the current nesting level
    istart = 0
    iend = 0
    call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))

    if (.not. bmatch) then
      ! Try to match the command
      call cmdprs_commandmatch (scommandtrimmed,spatternINT,bmatch)

      if (bmatch) then      
        ! Create INT variable.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,scommandtrimmed(istart:iend),inestlevel,ssectionname)
        call collct_setvalue_int (rcmdStatus%rcollection, scommandtrimmed(istart:iend), &
            0, .true., ssectionname=ssectionname) 
      end if
    end if
    
    if (.not. bmatch) then
      call cmdprs_commandmatch (scommandtrimmed,spatternDOUBLE,bmatch)
    
      if (bmatch) then
        ! Create DOUBLE variable.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,scommandtrimmed(istart:iend),inestlevel,ssectionname)
        call collct_setvalue_real (rcmdStatus%rcollection, scommandtrimmed(istart:iend), &
            0.0_DP, .true., ssectionname=ssectionname) 
      end if
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommandtrimmed,spatternSTRING,bmatch)
    
      if (bmatch) then
        ! Create STRING variable.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,scommandtrimmed(istart:iend),inestlevel,ssectionname)
        call collct_setvalue_string (rcmdStatus%rcollection, scommandtrimmed(istart:iend), &
            "", .true., ssectionname=ssectionname) 
      end if
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommandtrimmed,spatternINT2,bmatch)
    
      if (bmatch) then
        ! Create INT variable and assign.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,scommandtrimmed(istart:iend),inestlevel,ssectionname)
        call collct_setvalue_int (rcmdStatus%rcollection, scommandtrimmed(istart:iend), &
            0, .true., ssectionname=ssectionname)
        istart2 = 0
        iend2 = 0
        call tpsym_evalExpression (scommandtrimmed(istart:),rcmdStatus,inestlevel,istart2,iend2,rvalue)
      end if
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommandtrimmed,spatternDOUBLE2,bmatch)
    
      if (bmatch) then
        ! Create DOUBLE variable and assign.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,scommandtrimmed(istart:iend),inestlevel,ssectionname)
        call collct_setvalue_real (rcmdStatus%rcollection, scommandtrimmed(istart:iend), &
            0.0_DP, .true., ssectionname=ssectionname) 
        istart2 = 0
        iend2 = 0
        call tpsym_evalExpression (scommandtrimmed(istart:),rcmdStatus,inestlevel,istart2,iend2,rvalue)
      end if
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommandtrimmed,spatternSTRING2,bmatch)
    
      if (bmatch) then
        ! Create STRING variable and assign.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_nexttoken (scommandtrimmed,istart,iend,len(scommandtrimmed))
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,scommandtrimmed(istart:iend),inestlevel,ssectionname)
        call collct_setvalue_string (rcmdStatus%rcollection, scommandtrimmed(istart:iend), &
            "", .true., ssectionname=ssectionname) 
        istart2 = 0
        iend2 = 0
        call tpsym_evalExpression (scommandtrimmed(istart:),rcmdStatus,inestlevel,istart2,iend2,rvalue)
      end if
    end if

    if (bworkflowAllowed .and. (rcmdStatus%ierror .eq. 0)) then
      
      ! 2nd chance. Control commands parse the original string, must not contain a ";"
      ! at the end.

      if (.not. bmatch) then
        
        call cmdprs_commandmatch (scommandtrimmed,spatternBEGIN,bmatch)
      
        if (bmatch) then
          ! Get the subblock and execute
          call cmdprs_getNextBlock (rcmdblock,iline,istart2,iend2)
          
          ! Create new nest level in the collection
          call collct_addsection (rcmdStatus%rcollection, trim(sys_siL(inestlevel+1,10)))
          
          ! Execute
          call cmdprs_parsecmdblock (rcmdStatus,rvalue,rcmdblock,inestlevel+1,istart2,iend2)

          ! Remove local symbols
          call collct_deletesection (rcmdStatus%rcollection, trim(sys_siL(inestlevel+1,10)))

          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if
          
      if (.not. bmatch) then

        call cmdprs_commandmatch (scommandtrimmed,spatternFOR,bmatch)
      
        if (bmatch) then
          ! For loop. Fetch the command block from the next lines.
          call cmdprs_getNextBlock (rcmdblock,iline+1,istart2,iend2)
          
          ! Process the command
          call cmdprs_doFor (rcmdStatus,rcmdblock,inestlevel,scommandtrimmed,iline,istart2,iend2)
          
          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if
          
      if (.not. bmatch) then

        call cmdprs_commandmatch (scommandtrimmed,spatternWHILE,bmatch)

        if (bmatch) then
          ! While loop. Fetch the command block from the next lines.
          call cmdprs_getNextBlock (rcmdblock,iline+1,istart2,iend2)

          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if

      if (.not. bmatch) then

        call cmdprs_commandmatch (scommandtrimmed,spatternDO,bmatch)
      
        if (bmatch) then
          ! Do-while loop. Fetch the command block from the next lines.
          call cmdprs_getNextBlock (rcmdblock,iline+1,istart2,iend2)

          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if

      if (.not. bmatch) then

        call cmdprs_commandmatch (scommandtrimmed,spatternIF,bmatch)
      
        if (bmatch) then
          ! If-command. Most complicated. Get one or two subblocks.
          call cmdprs_getNextBlock (rcmdblock,iline+1,istart2,iend2)
          
          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if

    end if

    if (.not. bmatch) then
      
      ! No build-in command. Process special command.
      ! These may miss a training ";"!
      call cmdprs_dospecialcommand (rcmdStatus,scommandtrimmed,ierror)
      bmatch = ierror .eq. 0
      if (ierror .gt. 1) then
        ! =0: ok, =1: not found. >2: real error, pass to caller
        rcmdStatus%ierror = ierror
        bmatch = .true.
      end if
      
    end if
    
    if (.not. bmatch) then
!      if (.not. btrailexists) then
!        call output_line ("Error. "";"" missing at the end!")
!        rvalue%ctype = STYPE_INVALID
!        rcmdStatus%ierror = 1
!        return
!      end if

      ! Try it as an immediate expression.
      ! Processes e.g. variable assignments etc...
      istart = 0
      iend = 0
      call tpsym_evalExpression (scommandtrimmed,rcmdStatus,inestlevel,istart,iend,rvalue)
      bmatch = rvalue%ctype .ne. STYPE_INVALID
    end if
    
    if (.not. bmatch) then
      call output_line ("Invalid command!")
    else if (rcmdStatus%ierror .ne. 0) then
      call output_line ("Invalid syntax!")
    end if
      
  end subroutine

  ! ***************************************************************************

  recursive subroutine cmdprs_parsecmdblock (rcmdStatus,&
      rreturnValue,rcmdblock,inestlevel,iblockstart,iblockend)
  
  !<description>
    ! Parses a command block.
  !</description>
  
  !<inputoutput>
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command block with commands.
    ! The block may start with "{"; in this case it must end with "}".
    ! the block must at least contain one line.
    type(t_commandBlock), intent(in) :: rcmdblock
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Start line in the command block. If not specified, the first line is assumed.
    integer, intent(in), optional :: iblockstart
    
    ! End line in the command block. If not specified, the last line is assumed.
    integer, intent(in), optional :: iblockend
    
    ! Structure encapsuling return values.
    ! Initialised with default values by "intent(out)".
    type(t_symbolValue), intent(out) :: rreturnvalue
  !</input>

    integer :: iline, istart, iend, iblock
    character(len=SYS_STRLEN*2) :: sline,stemp2
    logical :: bfirst
    
    istart = 1
    iend = rcmdblock%icount
    if (present(iblockstart)) istart = iblockstart
    if (present(iblockend)) iend = iblockend
    
    ! Parse all lines.
    iline = 0
    iblock = 0
    bfirst = .true.
    do
    
      ! Get the line and split + dequote it.
      call cmdprs_getNextLine (rcmdblock,iline,istart,iend,sline)
      
      ! Stop if there is no more line or we have to stop the processing.
      if ((iline .eq. 0) .or. &
          rcmdStatus%bterminate  .or. (rcmdStatus%ierror .ne. 0)) exit
    
      if (rcmdStatus%becho) then
        stemp2 = rcmdblock%p_Rcommands(iline)
        call output_line (trim(stemp2))
      end if

      ! Is this a block start and the beginning of the block?
      if ((sline .eq. "{") .and. bfirst) then
        ! Start of our block. Increase the block counter.
        iblock = iblock + 1
      else if (sline .eq. "}")then
        ! Decrease the block counter
        iblock = iblock - 1
      else
        ! Execute this as a command.
        call cmdprs_docommand (rcmdStatus,rcmdblock,inestlevel,&
            sline,iline,iblockstart,iblockend,.true.,rreturnvalue)
      end if

      bfirst = .false.
    
    end do
    
    ! Error if the block end cannot be found.
    if (iblock .gt. 0) then
      call output_line ("Cannot find end of block!")
      rcmdStatus%ierror = 1
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_parsestream (rcmdStatus,istream,inestlevel)
  
  !<description>
    ! Parses an input stream.
  !</description>
  
  !<inputoutput>
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! The input stream channel.
    integer, intent(in) :: istream

    ! Level of nesting
    integer, intent(in) :: inestlevel
  !</input>

    type(VARYING_STRING) :: sinput,sin
    type(t_commandBlock) :: rcmdBlock
    integer :: ierr
    type(t_symbolValue) :: rreturnvalue

    ! Initialise the command block
    call cmdprs_initcmdblock (rcmdBlock)

    do while (.not. rcmdStatus%bterminate)
    
      sin = ""
      sinput = ""
      
      call get(istream,sin,IOSTAT=ierr) ! read next line of file
      
      ! Trim. Warning, this allocates memory!!!
      sinput = trim(sin)
      
      if (ierr == -1 .or. ierr > 0) then
    
        if (rcmdblock%icount .eq. 0) then
          ! Finish.
          exit
        end if
        
        ! There's no more data. Start executing.
        call cmdprs_parsecmdblock (rcmdStatus,rreturnvalue,rcmdblock,inestlevel)
  
        ! Create a new command block
        call cmdprs_donecmdblock (rcmdBlock)
        call cmdprs_initcmdblock (rcmdBlock)
        
      else
        if (len(sinput) .ne. 0) then
          ! Ignore comments
          if (getchar (sinput,1) .ne. "#") then
      
            ! Append the line to the command block
            call cmdprs_addCommand (rcmdBlock,sinput)
              
          end if
        end if
      end if
    end do

    ! Release memory    
    call cmdprs_donecmdblock (rcmdBlock)
    sin = ""
    sinput = ""
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_parseterminal (rcmdStatus)
  
  !<description>
    ! Parses an the terminal.
  !</description>

  !<inputoutput>
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

    type(VARYING_STRING) :: sinput
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
    character(len=SYS_STRLEN) :: sstr
    type(t_commandBlock) :: rcmdBlock
    type(t_symbolValue) :: rreturnvalue

    do while (.not. rcmdStatus%bterminate)
    
      ! Reset the error marker.
      rcmdStatus%ierror = 0
      
      ! Read tghe next line from the termninal
      call output_line("> ",bnolinebreak=.true.,bnotrim=.true.)
      read (IN_TERMINAL,'(A)') sstr
      
      sinput = trim(adjustl(sstr))
      
      if (rcmdStatus%becho) then
        call output_line (CHAR(sinput))
      end if
      
      if (len(sinput) .gt. 1) then
      
        ! Ignore comments
        if (getchar (sinput,1) .ne. "#") then

          ! Form a command block containing one command.
          call cmdprs_initcmdblock (rcmdBlock)
          call cmdprs_addCommand (rcmdBlock,sinput)
          
          ! Execute the command block
          call cmdprs_parsecmdblock (rcmdStatus,rreturnvalue,rcmdblock,0)
            
          ! Release the command block
          call cmdprs_donecmdblock (rcmdBlock)
          
        end if
      
      end if
      
      call output_lbrk()
     
    end do
    
    ! Release memory.
    sinput = ""
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine cmdprs_dospecialcommand (rcmdStatus,scommand,ierror)
  
  !<description>
    ! Executes a command.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Command string.
    character(len=*), intent(in) :: scommand
  !</inputoutput>
  
  !<output>
    ! Error flag.
    ! =0: no error.
    ! =1: Command not found
    integer, intent(out) :: ierror
  !</output>

    character(len=SYS_NAMELEN) :: sargument
    character(len=SYS_STRLEN) :: snewcommand,stemp
    integer :: istart,iend,ilength
    type(t_symbolValue) :: rvalue
    
    ierror = 0
    
    ! Get the command
    istart = 0
    iend = 0
    ilength = len_trim(scommand)
    call cmdprs_nexttoken (scommand,istart,iend,ilength)
    
    ! There must be no "("!
    call cmdprs_nexttoken (scommand,istart,iend,ilength)
    if (istart .gt. 0) then
      if (scommand(istart:iend) .eq. "(") then
        ierror = 1
        return
      end if
    end if
    
    call cmdprs_nexttoken (scommand,istart,iend,ilength,.true.)
    
    ! For some commands, formulate the correct syntax and evaluate the
    ! expression behind it.
    snewcommand = ""

    if (scommand(istart:iend) .eq. "exit") then
      snewcommand = "halt();"
    
    else if (scommand(istart:iend) .eq. "print") then
      call cmdprs_nexttoken (scommand,istart,iend,ilength)
      if (istart .ne. 0) then
        sargument = scommand(istart:iend)
        snewcommand = "printf(""%s"","//trim(sargument)//");"
      else
        snewcommand = "printf();"
      end if

    else if (scommand(istart:iend) .eq. "help") then
      call cmdprs_nexttoken (scommand,istart,iend,ilength)
      if (istart .ne. 0) then
        call cmdprs_dequoteStd (scommand(istart:iend),sargument)
        snewcommand = "help("""//trim(sargument)//""");"
      else
        snewcommand = "help();"
      end if

    else if (scommand(istart:iend) .eq. "run") then
      call cmdprs_nexttoken (scommand,istart,iend,ilength)
      if (istart .ne. 0) then
        call cmdprs_dequoteStd (scommand(istart:iend),sargument)
        snewcommand = "run("""//trim(sargument)//""");"
      else
        snewcommand = "run();"
      end if
    end if
    
    ! If a command was created, invoke it.
    if (snewcommand .ne. "") then
      call cmdprs_complexSplit (snewcommand,stemp,1,.false.)
      istart = 0
      iend = 0
      call tpsym_evalExpression (stemp,rcmdStatus,0,istart,iend,rvalue)
      return
    end if    

    ! Not found
    ierror = 1
  
  end subroutine
  
  ! ***************************************************************************
  ! Parsing of symbols
  ! ***************************************************************************
  
  ! ***************************************************************************

  !<subroutine>
  
  recursive subroutine tpsym_evalExpression (sstring,rcmdStatus,inestlevel,istart,iend,rvalue,&
      spreviousOperator)
  
  !<description>
    ! Parses a string of symbols and creates a value from it.
  !</description>

  !<input>
    ! String containing symbols. The string must have been split with cmdprs_complexSplit
    ! and must completely be trimmed!
    character(len=*), intent(in) :: sstring
  !</input>

  !<inputoutput>
    ! IN: Start/End of the first token of the expression or =0,
    ! if to start with the very first token.
    ! OUT: Start/End of the next token after the expression.
    ! =0 if the expression is completed.
    integer, intent(inout) :: istart,iend
  !</inputoutput>

  !<input>
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Level of nesting
    integer, intent(in) :: inestlevel
    
    ! OPTIONAL: The operator preceding this expression. Must be specified
    ! in sequences of expressions. If omitted, the value of the expression
    ! is just calculated and returned, ignoring any operator precedence.
    character(len=*), intent(in), optional :: spreviousOperator
  !</input>
  
  !<input
  
  !<output>
    ! The value of the expression.
    type(t_symbolValue), intent(out) :: rvalue
  !</output>
  
  !</subroutine>
  
    ! local variables
    integer :: icurropstart,icurropend,isymbolstart,isymbolend
    integer :: inextsymstart, inextsymend
    type(t_symbolValue) :: rvaluetemp,rvaluemultiplier
    integer :: ilength, ibracketcount
    logical :: berror
    type(t_symbolValue), dimension(:), pointer :: p_Rvalues
    
    ilength = len_trim(sstring)
    
    ! Start with the first token?
    if (istart .eq. 0) then
      call cmdprs_nexttoken (sstring,istart,iend,ilength)
    end if

    ! Current token.
    isymbolstart = istart
    isymbolend = iend
    
    ! How does the expression start?
    if (sstring(istart:iend) .eq. "(") then
      ! Evaluate the value in the brackets.
      call tpsym_evalExpressionBrackets (sstring,rcmdStatus,inestlevel,istart,iend,rvalue)
      
      ! istart/iend now points to the next token after the closing bracket.
    else if (sstring(isymbolstart:isymbolend) .eq. "-") then
    
      ! This negates the whole value. Evaluate the rest, then multiply by -1.
      ! Cancel if there is a previous operator -- we do not allow something like "* -1".
      if (.not. present(spreviousOperator)) then
        ! Ignore the "-".
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
      
        ! Evaluate the next stuff
        call tpsym_evalExpression (sstring,rcmdStatus,inestlevel,istart,iend,rvalue)
            
        ! Multiply by -1.
        rvaluetemp%ctype = STYPE_INTEGER
        rvaluetemp%ivalue = -1
        rvalue = rvalue * rvaluetemp
      
      end if

    else
      ! Check the following token. If it's a "(", we have a function call!
      inextsymstart = istart
      inextsymend = iend
      call cmdprs_nexttoken (sstring,inextsymstart,inextsymend,ilength)
      if (inextsymstart .eq. 0) then
        ! Variable or symbol. Parse the symbol, create a value.
        call tpsym_parseSymbol (sstring(isymbolstart:isymbolend),rcmdStatus%rcollection,inestlevel,rvalue)
        
        ! Go to the next token
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
        
      else if (sstring(inextsymstart:inextsymend) .eq. "(") then
      
        istart = inextsymstart
        iend = inextsymend
        
        ! That's a function. Get the arguments.
        call tpsym_getArguments (sstring,rcmdStatus,inestlevel,istart,iend,p_Rvalues,berror)
        
        if (.not. berror) then
          ! Evaluate the function behind it.
          if (associated(p_Rvalues)) then
            call tpsym_evalFunction (sstring(isymbolstart:isymbolend),rcmdStatus,inestlevel,&
                rvalue,p_Rvalues)

            ! Release memory.
            deallocate(p_Rvalues)
          else
            call tpsym_evalFunction (sstring(isymbolstart:isymbolend),rcmdStatus,inestlevel,&
                rvalue)
          end if
          
        end if
        
        ! istart/iend points to the next token.
      else
        ! Variable or symbol. Parse the symbol, create a value.
        call tpsym_parseSymbol (sstring(isymbolstart:isymbolend),&
            rcmdStatus%rcollection,inestlevel,rvalue)
        
        ! Go to the next token
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
        
      end if
    end if
    
    ! Return if invalid or variable which cannot be evaluated.
    if (rvalue%ctype .eq. STYPE_INVALID) then
      ! Return an error.
      rcmdStatus%ierror = 1
      rvalue%ctype = STYPE_INTEGER
      return
    end if
    if (rvalue%ctype .eq. STYPE_VAR) return
    
    ! Repeat to evaluate all following operators.
    do while (istart .ne. 0)
      ! We have a value. Get the next operator.
      icurropstart = istart
      icurropend = iend
      
      ! What is more important, the left or the rigt operator?
      !
      ! Left more important -> Stop here, return, value will be used.
      ! Right more important -> Recursively evaluate the expression on the right
      !    and evaluate the value of the operator.
      if (present(spreviousOperator)) then
        if (icurropstart .eq. 0) return
        if (getOperatorWeight(spreviousOperator) .gt. &
            getOperatorWeight(sstring(icurropstart:icurropend))) return
      end if
      
      ! Check for immediate operators.
      if (sstring(icurropstart:icurropend) .eq. ",") then
        ! Cancel, we are done.
        return
        
      else if (sstring(icurropstart:icurropend) .eq. ";") then
        ! Cancel, we are done.
        return
        
      else if (sstring(icurropstart:icurropend) .eq. ")") then
        ! Cancel, we are done.
        return
        
      else if (sstring(icurropstart:icurropend) .eq. "++") then
      
        ! Increase variable, return old value.
        rvaluetemp = 1
        rvalue = rvalue + rvaluetemp
        call tpsym_saveSymbol (rvalue)
        rvalue = rvalue - rvaluetemp
        
        ! Next token
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
        if (istart .eq. 0) return
        
      elseif (sstring(icurropstart:icurropend) .eq. "--") then
      
        ! Increase variable, return old value.
        rvaluetemp = 1
        rvalue = rvalue - rvaluetemp
        call tpsym_saveSymbol (rvalue)
        rvalue = rvalue + rvaluetemp

        ! Next token
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
        if (istart .eq. 0) return
        
      else      

        ! Right operator more important. Evaluate the right expression first and
        ! calculate a new value with the operator.
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
        if (istart .eq. 0) return
        
        call tpsym_evalExpression (sstring,rcmdStatus,inestlevel,istart,iend,rvaluetemp,&
            spreviousOperator = sstring(icurropstart:icurropend))
        if (rvaluetemp%ctype .eq. STYPE_INVALID) return
        
        ! What is the operator?
        if (sstring(icurropstart:icurropend) .eq. "+") then
          rvalue = rvalue + rvaluetemp
        else if (sstring(icurropstart:icurropend) .eq. "-") then
          rvalue = rvalue - rvaluetemp
        else if (sstring(icurropstart:icurropend) .eq. "*") then
          rvalue = rvalue * rvaluetemp
        else if (sstring(icurropstart:icurropend) .eq. "/") then
          rvalue = rvalue / rvaluetemp
        else if (sstring(icurropstart:icurropend) .eq. "%") then
          rvalue = mod(rvalue, rvaluetemp)
        else if (sstring(icurropstart:icurropend) .eq. "==") then
          ! Compare. Result is an int without any connection to a variable.
          if (rvalue == rvaluetemp) then
            call tpsym_undefine(rvalue)
            rvalue = 1
          else
            rvalue = 0
          end if
        else if (sstring(icurropstart:icurropend) .eq. "!=") then
          ! Compare. Result is an int without any connection to a variable.
          if (rvalue /= rvaluetemp) then
            call tpsym_undefine(rvalue)
            rvalue = 1
          else
            call tpsym_undefine(rvalue)
            rvalue = 0
          end if
        else if (sstring(icurropstart:icurropend) .eq. "<") then
          ! Compare. Result is an int without any connection to a variable.
          if (rvalue < rvaluetemp) then
            call tpsym_undefine(rvalue)
            rvalue = 1
          else
            call tpsym_undefine(rvalue)
            rvalue = 0
          end if
        else if (sstring(icurropstart:icurropend) .eq. ">") then
          ! Compare. Result is an int without any connection to a variable.
          if (rvalue > rvaluetemp) then
            call tpsym_undefine(rvalue)
            rvalue = 1
          else
            call tpsym_undefine(rvalue)
            rvalue = 0
          end if
        else if (sstring(icurropstart:icurropend) .eq. "<=") then
          ! Compare. Result is an int without any connection to a variable.
          if (rvalue <= rvaluetemp) then
            call tpsym_undefine(rvalue)
            rvalue = 1
          else
            call tpsym_undefine(rvalue)
            rvalue = 0
          end if
        else if (sstring(icurropstart:icurropend) .eq. ">=") then
          ! Compare. Result is an int without any connection to a variable.
          if (rvalue >= rvaluetemp) then
            call tpsym_undefine(rvalue)
            rvalue = 1
          else
            call tpsym_undefine(rvalue)
            rvalue = 0
          end if
        else if (sstring(icurropstart:icurropend) .eq. "=") then
          ! Assign + save.
          rvalue = rvaluetemp
          call tpsym_saveSymbol (rvalue)
        end if
        
      end if

    end do
    
  contains
  
    ! --------------------------------------------
    integer function getOperatorWeight (soperator)
    
    ! Returns the weight of an operator.
    ! Use the C operator precedence:
    !
    !    http://www.difranco.net/cop2220/op-prec.htm

    ! The operator
    character(len=*), intent(in) :: soperator
    
      getOperatorWeight = 0
      
      if (soperator .eq. "||") then
        getOperatorWeight = 3
      else if (soperator .eq. "&&") then
        getOperatorWeight = 4
      else if (soperator .eq. "|") then
        getOperatorWeight = 5
      else if (soperator .eq. "^") then
        getOperatorWeight = 6
      else if (soperator .eq. "&") then
        getOperatorWeight = 7
      else if ((soperator .eq. "==") .or. (soperator .eq. "!=")) then
        getOperatorWeight = 8
      else if ((soperator .eq. "<") .or. (soperator .eq. "<=") .or. &
               (soperator .eq. ">") .or. (soperator .eq. ">=")) then
        getOperatorWeight = 9
      else if ((soperator .eq. "<<") .or. (soperator .eq. ">>")) then
        getOperatorWeight = 10
      else if ((soperator .eq. "+") .or. (soperator .eq. "-")) then
        getOperatorWeight = 11
      else if ((soperator .eq. "*") .or. (soperator .eq. "/") .or. &
               (soperator .eq. "%")) then
        getOperatorWeight = 12
      else if ((soperator .eq. "++") .or. (soperator .eq. "--") .or. &
               (soperator .eq. "!") .or. (soperator .eq. "~")) then
        getOperatorWeight = 13
      end if
    
    end function
  
  end subroutine

  ! ***************************************************************************

  !<subroutine>
  
  recursive subroutine tpsym_evalExpressionBrackets (sstring,rcmdStatus,inestlevel,istart,iend,rvalue)
  
  !<description>
    ! Parses an expression in brackets.
  !</description>

  !<input>
    ! String containing symbols. The string must have been split with cmdprs_complexSplit
    ! and must completely be trimmed!
    character(len=*), intent(in) :: sstring
  
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Level of nesting
    integer, intent(in) :: inestlevel
  !</input>
  
  !<inputoutput>
    ! IN: Start/End of the opening bracket.
    ! OUT: Start/End of the closing bracket.
    integer, intent(inout) :: istart,iend
  !</inputoutput>
  
  !<output>
    ! The value of the expression.
    type(t_symbolValue), intent(out) :: rvalue
  !</output>
  
  !</subroutine>
  
    ! local variables
    integer :: istart2,iend2,iend3
    integer :: ibracket,ilength
    
    ! Get basic information
    ilength = len(sstring)
    
    ! Expression in brackets. Most important.
    !
    ! Search closing bracket -- or inner expression(s) to be evaluated.
    istart2 = istart
    iend2 = iend
    iend3 = iend
    ibracket = 1
    call cmdprs_nexttoken (sstring,istart2,iend2,ilength)
    do while ((ibracket .gt. 0) .and. (istart2 .ne. 0))
      if (sstring(istart2:istart2) .eq. "(") then
        ibracket = ibracket + 1
      else if (sstring(istart2:istart2) .eq. ")") then
        ibracket = ibracket - 1
      end if
      
      ! Ignore other tokens.
      call cmdprs_nexttoken (sstring,istart2,iend2,ilength)
    end do
    
    ! Return if ibracket <> 0, expression invalid.
    if (ibracket .gt. 0) return
    
    ! Evaluate recursively.
    istart = 0
    iend = 0
    call tpsym_evalExpression (sstring(iend3+2:istart2-4),rcmdStatus,inestlevel,&
        istart,iend,rvalue)
    
    ! Hop behind the bracket.
    istart = istart2
    iend = iend2

  end subroutine
  
  ! ***************************************************************************

  !<subroutine>
  
  recursive subroutine tpsym_getArguments (ssource,rcmdStatus,inestlevel,istart,iend,p_Rvalues,berror)
  
  !<description>
    ! Evaluates a list of arguments in brackets.
  !</description>

  !<input>
    ! String containing symbols. The string must have been split with cmdprs_complexSplit
    ! and must completely be trimmed!
    character(len=*), intent(in) :: ssource
  
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Level of nesting
    integer, intent(in) :: inestlevel
  !</input>
  
  !<inputoutput>
    ! IN: Start/End of the function name.
    ! OUT: Start/End of the token behind the closing bracket.
    integer, intent(inout) :: istart,iend
  !</inputoutput>
  
  !<output>
    ! Pointer to list of values or NULL if the list is empty.
    type(t_symbolValue), dimension(:), pointer :: p_Rvalues
    
    ! Set to TRUE in case of a format error
    logical, intent(out) :: berror
  !</output>
  
  !</subroutine>
  
    ! local variables
    integer :: istart2,iend2
    integer :: iargs,ilength,i,ibracket
    
    berror = .false.

    ! ignore the bracket.
    ilength = len_trim(ssource)
    call cmdprs_nexttoken (ssource,istart,iend,ilength)
    
    ! Count the number of "," before teh closing bracket.
    istart2 = istart
    iend2 = iend
    iargs = 0
    if (ssource(istart:iend) .ne. ")") then
      iargs = 1
      ibracket = 1
      do
        if (istart .eq. 0) return
        if (ssource(istart:iend) .eq. "(") then
          ibracket = ibracket+1
        else if (ssource(istart:iend) .eq. ")") then
          ibracket = ibracket-1
          ! Check if we are done.
          if (ibracket .eq. 0) exit 
        else if ((ssource(istart:iend) .eq. ",") .and. (ibracket .eq. 1)) then  
          ! Our argument
          iargs = iargs + 1
        end if
        
        ! Next token
        call cmdprs_nexttoken (ssource,istart,iend,ilength)
      end do
      if (ibracket .ne. 0) then
        ! Error
        berror = .true.
        return
      end if
    end if
    
    ! Allocate memory for the subexpressions.
    if (iargs .gt. 0) then
      allocate(p_Rvalues(iargs))
      
      ! Get the subexpressions
      istart = istart2
      iend = iend2
      do i=1,iargs
      
        ! Take a look at the current and next token. If we have the form
        ! "identifier = ...", this is the specification of an optional
        ! variable.
        istart2 = istart
        iend2 = iend
        call cmdprs_nexttoken (ssource,istart2,iend2,ilength)
        if (istart2 .ne. 0) then
          if (ssource(istart2:iend2) .eq. "=") then
            ! Skip the name as well as the "=" sign.
            call cmdprs_nexttoken (ssource,istart,iend,ilength)
            call cmdprs_nexttoken (ssource,istart,iend,ilength)
            
            ! Then evaluate...
          end if
        end if
      
        ! Evaluate the next expression up to the closing bracket.
        call tpsym_evalExpression (ssource,&
            rcmdStatus,inestlevel,istart,iend,p_Rvalues(i))            
            
        if (istart2 .ne. 0) then            
          if (ssource(istart2:iend2) .eq. "=") then
            ! Optional argument. Get its name as variable tag.
            ! For that purpose, go one token back from the "=" sign.
            call cmdprs_nexttoken (ssource,istart2,iend2,ilength,.true.)
            p_Rvalues(i)%svartag = ssource(istart2:iend2)
          end if
        end if
        
        ! istart/iend point to the "," or ")".
        ! Ignore.
        call cmdprs_nexttoken (ssource,istart,iend,ilength)
      
      end do
      
    else
      nullify(p_Rvalues)
    end if

  end subroutine

  ! ***************************************************************************
  ! Work flow command handlers
  ! ***************************************************************************

  ! ***************************************************************************

  !<subroutine>

  subroutine cmdprs_doFor (rcmdStatus,rcmdblock,inestlevel,scommand,iline,&
      iblockstart,iblockend)
  
  !<description>
    ! Work flow command: FOR loop
  !</description>

  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>
  
  !<input>
    ! Current command block containing the program
    type(t_commandBlock), intent(in) :: rcmdBlock

    ! Command string. Splitted.
    character(len=*), intent(in) :: scommand
    
    ! Number of the line
    integer, intent(inout) :: iline
    
    ! Start of the commands depending on the command in rcmdBlock.
    integer, intent(in) :: iblockstart

    ! End of the commands depending on the command in rcmdBlock.
    integer, intent(in) :: iblockend
    
    ! Level of nesting
    integer, intent(in) :: inestlevel
  !</input>
  
  !</subroutine>

    ! local variables
    integer :: isemicolon1,isemicolon2
    integer :: istart,iend,ilength,isubcommandline
    integer :: istartarg, iendarg
    logical :: berror
    logical :: bfirst
    type(t_symbolValue) :: rreturnvalue

    ! scommand has a special stucture:
    !   FOR ( ... )
    ! The entries in the braces have the form
    !    ... ; ... ; ...
    ! with "..." specifying an expression.
    
    ! Search for the two semicola
    istart = 0
    ilength = len_trim(scommand)
    isemicolon1 = 0
    isemicolon2 = 0
    berror = .false.
    call cmdprs_nexttoken (scommand,istart,iend,ilength)
    do while (istart .ne. 0) 
      ! Find the semicola
      if (scommand(istart:iend) .eq. ";") then
        if (isemicolon1 .eq. 0) then
          isemicolon1 = istart
        else if (isemicolon2 .eq. 0) then
          isemicolon2 = istart
        else
          berror = .true.
        end if

      ! Find the start of the arguments (behind the "(") and
      ! their ends (before the ")").
      else if (scommand(istart:iend) .eq. "(") then
        istartarg = iend+2
      else if (scommand(istart:iend) .eq. ")") then  
        iendarg = istart-2
      end if

      call cmdprs_nexttoken (scommand,istart,iend,ilength)
    end do
    if (isemicolon2 .eq. 0) berror = .true.
    
    if (.not. berror) then
      ! Loop initialisation
      isubcommandline = iline+1
      call cmdprs_docommand (rcmdStatus,rcmdblock,inestlevel,&
          scommand(istartarg:isemicolon1-2),isubcommandline,&
          iblockstart,iblockend,.false.,rreturnvalue)
    
      ! Execute the FOR loop in an endless DO loop.
      bfirst = .true.
      do while (.not. rcmdStatus%bterminate .and. (rcmdStatus%ierror .eq. 0)) 
      
        if (.not. bfirst) then
          ! Check termination criterion
          isubcommandline = iline+1
          call cmdprs_docommand (rcmdStatus,rcmdblock,inestlevel,&
              scommand(isemicolon1+2:isemicolon2-2),isubcommandline,&
              iblockstart,iblockend,.false.,rreturnvalue)
              
          ! Stop if FALSE is returned.
          if (rreturnvalue%ivalue .eq. 0) exit
          
        end if
        
        if (.not. rcmdStatus%bterminate .and. (rcmdStatus%ierror .eq. 0)) then
        
          ! The actual loop commands.
          call cmdprs_parsecmdblock (rcmdStatus,rreturnvalue,rcmdblock,inestlevel,&
              iblockstart,iblockend)

        end if

        if (.not. rcmdStatus%bterminate .and. (rcmdStatus%ierror .eq. 0)) then
          ! Final command block of the FOR loop.
          isubcommandline = iline+1
          call cmdprs_docommand (rcmdStatus,rcmdblock,inestlevel,&
              scommand(isemicolon2+2:iendarg),isubcommandline,&
              iblockstart,iblockend,.false.,rreturnvalue)
        end if
        
        bfirst = .false.
      
      end do
    end if
    
    if (berror) then
      call output_line ("Wrong syntax!")
    end if

  end subroutine

  ! ***************************************************************************
  ! Extended command handlers
  ! ***************************************************************************

  !<subroutine>
  
  recursive subroutine tpsym_evalFunction (sstring,rcmdStatus,inestlevel,rvalue,Rvalues)
  
  !<description>
    ! Evaluates a function.
  !</description>

  !<input>
    ! Name of the function
    character(len=*), intent(in) :: sstring
  
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! OPTIONAL: list of values or NULL if the list is empty.
    type(t_symbolValue), dimension(:), optional :: Rvalues
  !</input>
  
  !<output>
    ! Return value of the function.
    type(t_symbolValue), intent(out) :: rvalue
  !</output>
  
  !</subroutine>
  
    character (len=SYS_STRLEN) :: stemp
    logical :: bunknown
    integer :: iunit
    logical :: bexists

    rvalue%ctype = STYPE_INVALID
    
    bunknown = .false.
    
    if (sstring .eq. "run") then
    
      if (present(Rvalues)) then
        if (size(Rvalues) .eq. 1) then
          if (Rvalues(1)%ctype .eq. STYPE_STRING) then
          
            ! Open the file and read.
            call cmdprs_dequoteStd(Rvalues(1)%svalue,stemp)
            inquire(file=trim(stemp), exist=bexists)
            if (.not. bexists) then
              call output_line ("File not found!")
            else
              
              ! Open
              call io_openFileForReading(trim(stemp), iunit, .true.)

              ! Interpret the stream
              call cmdprs_parsestream (rcmdStatus,iunit,0)
              
              close (iunit)
            end if
            
            ! Worked.
            rvalue%ctype = STYPE_INTEGER
            return
          end if
        end if
        
      end if

      ! If we come to here, there is something wrong.
      call output_line ("Invalid arguments!")
      return

    else 
    
      ! Try to evaluate a FEAT command
      call fcmd_evalCommand (sstring,rcmdStatus,inestlevel,rvalue,bunknown,Rvalues)

    end if
    
    if (bunknown) then
      call output_line ("Command not found!")
      return
    end if

  end subroutine
  
end module
