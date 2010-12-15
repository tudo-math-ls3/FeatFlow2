module commandparser

  use ISO_VARYING_STRING
  
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

  implicit none
  
  private
  
  ! Type that encapsules the current status.
  type t_commandstatus
    
    ! Set to TRUE to terminate
    logical :: bterminate = .false.
    
    ! Error flag. =0: no error.
    integer :: ierror = 0
    
    ! Echo of the command line.
    logical :: becho = .false.
    
    ! Global collection object with all variables
    type(t_collection) :: rcollection
    
  end type
  
  ! Block of commands.
  type t_commandBlock
  
    ! Number of lines.
    integer :: icount = 0
  
    ! List of lines.
    type(VARYING_STRING), dimension(:), pointer :: p_Rcommands
  end type
  
  public :: t_commandstatus
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

  !************************************************************************

  elemental function cmdprs_dequote (ssource)
  !<description>
    ! Remove all quotation marks from a string.
  !</description>
  
  !<input>
    ! A source string.
    type(VARYING_STRING), intent(in) :: ssource
  !</input>
  
  !<output>
    ! Destination string. Non-escaped quotation marks are removed.
    type(VARYING_STRING) :: cmdprs_dequote
  !</output>
  
    type(VARYING_STRING) :: stemp,sdest
    integer :: i, idest
    character :: squotechar, scurrentchar
    logical :: bescape
    
    ! Copy the string, remove leading/trailing spaces.
    stemp = trim(ssource)
    
    ! Initialise destination
    sdest = stemp
    
    ! Loop through the string. Copy characters.
    ! Replace quotation marks.
    i = 1
    idest = 0
    bescape = .false.
    squotechar = ' '
    do 
      if (i .gt. len(stemp)) exit
      
      ! Get the current character.
      scurrentchar = getchar (stemp,i);
      
      if (bescape) then
      
        ! Switch off escape mode.
        bescape = .false.
        
      else
        
        if (scurrentchar .eq. "\") then
        
          ! take the next character as it is.
          bescape = .true.
          
          ! Ignore the escape char.
          i = i+1
          cycle
        
        else if ((scurrentchar .eq. "'") .or. (scurrentchar .eq. """")) then
          
          ! Start quoting or stop it. Only the current quote char
          ! triggers off the quoting.
          if (squotechar .eq. scurrentchar) then
            squotechar = " "
          else
            squotechar = scurrentchar
          end if
          
          ! Ignore the quotation char.
          i = i+1
          cycle
          
        end if
      end if

      ! Copy the character and go to the next char.
      idest = idest+1
      call setchar (sdest,idest,scurrentchar)

      i = i+1
    end do
    
    ! From the destination, take the actual substring.
    cmdprs_dequote = extract (sdest,1,idest)
    
    ! Release memory.
    sdest = ""
    stemp = ""
  
  end function 

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

  subroutine cmdprs_complexSplit (ssource,sdest,icharset,bdequote)
  
  !<description>
    ! Subduvides a line into an array of substrings.
    ! Quotation marks around substrings are removed.
    ! Escaped character are changed to standard characters.
    ! The different words are separated by "\0" characters.         
    ! the string is automatically trimmed.
  !</description>
  
  !<input>
    ! A source string.
    character(len=*), intent(in) :: ssource
    
    ! Type of character set to use for splitting.
    ! =0: split on word boundaries
    ! =1: split on expression boundaries (expression characters
    !     like brackets are separate words).
    ! =2: split on parameter boundaries (command line version)
    ! =3: spöit on parameter boundaries (sourcecode version)
    integer :: icharset

    ! Determins whether to de-quote substrings or not.
    logical, intent(in) :: bdequote
  !</input>

  !<output>
    ! The output string
    character(len=*), intent(out) :: sdest
  !</output>
  
    integer :: ipos,idest,imaxlen
    integer :: ccurrchargroup
    character :: ccurrentquote
    logical :: btriggered
    integer, parameter :: CHARGROUP_WHITESPACE = 1
    integer, parameter :: CHARGROUP_ONECHAREXP = 2
    integer, parameter :: CHARGROUP_TWOCHAREXP = 3
    integer, parameter :: CHARGROUP_INWORD     = 4
    integer, parameter :: CHARGROUP_INQUOTE    = 5
    
    ! The following chargroups are defined:
    ! Chargroup 1: Whitespaces. Ignored.
    ! Chargroup 2: One-character words. A character from this group
    !              triggers this character to be a word and the next character
    !              to be also a new word.
    ! Chargroup 3: Two-character words. This group contains actually two character
    !              groups. If the character from the first subgroup is found,
    !              The character from the 2nd subgroup must match. In this case,
    !              both characters together form a 2-character word.
    !              Examples: Expressions like "&&", "!=" or "==".
    ! Chargroup 4: Word characters. A word character behind a character from a
    !              lower chargroup triggers a new word. A "\" character triggers
    !              an escaped character which is transferred to the destination
    !              without change.
    ! Chargroup 5: Quotation marks with de-escape.
    !              1st occurence triggers start, 2nd occurrence the end.
    !              Escaped characters in this set are de-escaped.
    !              All characters in quotation marks are associated to
    !              this chargroup.
    ! Chargroup 6: Comment characters.
    !              These characters immediately stop the parsing.
    !
    ! The following character sets are defined:
    ! icharset=0: Splitting on word boundaries
    !    Group 1: " "
    !    Group 2: ""
    !    Group 3: "" / ""
    !    Group 4: Everything which is not in the other groups
    !    Group 5: "'""
    !    Group 6: "#"
    !
    ! icharset=1: Splitting on expression boundaries.
    !    Group 1: " "
    !    Group 2: "!+-*/%()~;{}^<>="
    !    Group 3: subgroup 1: "&|=!<>-+<>"
    !             subgroup 2: "&|==<>-+=="
    !    Group 4: Everything which is not in the other groups
    !    Group 5: "'""
    !    Group 6: "#"
    !
    ! icharset=2: Splitting on parameters boundaries (command line parameters)
    !    Group 1: " "
    !    Group 2: "="
    !    Group 3: "" / ""
    !    Group 4: Everything which is not in the other groups
    !    Group 5: "'""
    !    Group 6: "#"
    !
    ! icharset=3: Splitting on parameters boundaries (sourcecode parameters)
    !    Group 1: " "
    !    Group 2: "(),="
    !    Group 3: "" / ""
    !    Group 4: Everything which is not in the other groups
    !    Group 5: "'""
    !    Group 6: "#"
    !

    ! Clear the destination.
    sdest = ""

    ! Loop through the chars.
    ipos = 1
    idest = 0
    imaxlen = len_trim(ssource)
    
    ccurrchargroup = CHARGROUP_WHITESPACE
    btriggered = .false.
    
    do 
    
      ! Stop is wew are behind the last character.
      if (ipos .gt. imaxlen) exit
    
      ! Behaviour depends on current character group and set.
      select case (icharset)
      case (0)
        call trigger_chargroup5 (ccurrchargroup,"'""",btriggered,ssource,ipos,imaxlen,&
            sdest,idest,ccurrentquote,bdequote)
            
        if (.not. btriggered) then
          call trigger_chargroup1 (ccurrchargroup," ",btriggered,ssource,ipos)
        end if
        
      case (1)
        call trigger_chargroup5 (ccurrchargroup,"'""",btriggered,ssource,ipos,imaxlen,&
            sdest,idest,ccurrentquote,bdequote)
            
        if (.not. btriggered) then
          call trigger_chargroup3 (ccurrchargroup,"&|=!<>-+<>",&
                                                  "&|==<>-+==",&
              btriggered,ssource,ipos,imaxlen,sdest,idest)
        end if
        
        if (.not. btriggered) then
          call trigger_chargroup2 (ccurrchargroup,"!+-*/%()~;{}^<>=",btriggered,ssource,ipos,&
              sdest,idest)
        end if
        
        if (.not. btriggered) then
          call trigger_chargroup1 (ccurrchargroup," ",btriggered,ssource,ipos)
        end if
        
      case (2)
        call trigger_chargroup5 (ccurrchargroup,"'""",btriggered,ssource,ipos,imaxlen,&
            sdest,idest,ccurrentquote,bdequote)
            
        if (.not. btriggered) then
          call trigger_chargroup2 (ccurrchargroup,"=",btriggered,ssource,ipos,&
              sdest,idest)
        end if
            
        if (.not. btriggered) then
          call trigger_chargroup1 (ccurrchargroup," ",btriggered,ssource,ipos)
        end if
        
      case (3)
        call trigger_chargroup5 (ccurrchargroup,"'""",btriggered,ssource,ipos,imaxlen,&
            sdest,idest,ccurrentquote,bdequote)
            
        if (.not. btriggered) then
          call trigger_chargroup2 (ccurrchargroup,",=()",btriggered,ssource,ipos,&
              sdest,idest)
        end if
            
        if (.not. btriggered) then
          call trigger_chargroup1 (ccurrchargroup," ",btriggered,ssource,ipos)
        end if
        
      end select

      ! Trigger comment characters
      if (.not. btriggered) then
        call trigger_chargroup6 (ccurrchargroup,"#",btriggered,ssource,ipos,imaxlen)
      end if

      ! Remaining characters belong to group 4.
      if (.not. btriggered) then
        call process_chargroup4 (ccurrchargroup,ssource,ipos,imaxlen,sdest,idest)
      end if
    
    end do
    
  contains
  
    ! -------------------------------------------------------------------------
    subroutine trigger_chargroup1 (ccurrentgroup,schars,btriggered,ssource,isourcepos)
    
      ! Checks if the next token belongs to this character group.
      ! If that is the case, the destination and the pointers are modified according
      ! to the current state.
      
      ! Current character group
      integer, intent(inout) :: ccurrentgroup

      ! Characters in this group.
      character(len=*), intent(in) :: schars
      
      ! Returns if the group is triggered.
      logical, intent(out) :: btriggered
      
      ! Source string
      character(len=*), intent(in) :: ssource
      
      ! Source pointer
      integer, intent(inout) :: isourcepos
      
      ! local variables
      integer :: i
      
      i = index(schars,ssource(isourcepos:isourcepos))
      btriggered = i .ne. 0
      
      if (btriggered) then
        isourcepos = isourcepos+1

        ! We are now in whitespace-mode
        ccurrentgroup = CHARGROUP_WHITESPACE
      end if
    
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine trigger_chargroup2 (ccurrentgroup,schars,btriggered,ssource,isourcepos,&
        sdest,idestpos)
    
      ! Checks if the next token belongs to this character group.
      ! If that is the case, the destination and the pointers are modified according
      ! to the current state.
      
      ! Current character group
      integer, intent(inout) :: ccurrentgroup

      ! Characters in this group.
      character(len=*), intent(in) :: schars
      
      ! Returns if the group is triggered.
      logical, intent(out) :: btriggered
      
      ! Source string
      character(len=*), intent(in) :: ssource
      
      ! Source pointer
      integer, intent(inout) :: isourcepos
      
      ! Destination string
      character(len=*), intent(inout) :: sdest
      
      ! Destination pointer
      integer, intent(inout) :: idestpos
      
      ! local variables
      integer :: i
      
      i = index(schars,ssource(isourcepos:isourcepos))
      btriggered = i .ne. 0

      if (btriggered) then
        ! First character? probably insert a separator.
        if (idestpos .gt. 0) then
          idestpos = idestpos + 1;
          sdest(idestpos:idestpos) = char(0)
        end if

        ! Transfer the character
        idestpos = idestpos + 1
        sdest(idestpos:idestpos) = ssource(isourcepos:isourcepos)
        
        isourcepos = isourcepos+1
      
        ! We are now in one-character-expression mode
        ccurrentgroup = CHARGROUP_ONECHAREXP
      end if
    
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine trigger_chargroup3 (ccurrentgroup,schars1,schars2,btriggered,ssource,isourcepos,&
        imaxlen,sdest,idestpos)
    
      ! Checks if the next token belongs to this character group.
      ! If that is the case, the destination and the pointers are modified according
      ! to the current state.

      ! Current character group
      integer, intent(inout) :: ccurrentgroup
      
      ! Characters in this group. 1st subgroup.
      character(len=*), intent(in) :: schars1

      ! Characters in this group. 2nd subgroup.
      character(len=*), intent(in) :: schars2
      
      ! Returns if the group is triggered.
      logical, intent(out) :: btriggered
      
      ! Source string
      character(len=*), intent(in) :: ssource
      
      ! Source pointer
      integer, intent(inout) :: isourcepos
      
      ! Length of ssource
      integer, intent(inout) :: imaxlen
      
      ! Destination string
      character(len=*), intent(inout) :: sdest
      
      ! Destination pointer
      integer, intent(inout) :: idestpos
      
      ! local variables
      integer :: i,j
      
      i = 0
      do
        ! Find the character. This is a 2-character matching test, so if there
        ! is only 1 character left, we cannot match!
        ! Continue the search if the previous 2nd char did not match.
        j = index(schars1(i+1:),ssource(isourcepos:isourcepos))
        i = i+j
        btriggered = (j .ne. 0) .and. (isourcepos .lt. imaxlen)

        if (btriggered) then
          ! Check if the next character (if it exists) matches the character of
          ! the next subgroup.
          btriggered = schars2(i:i) .eq. ssource(isourcepos+1:isourcepos+1)
          
          if (btriggered) then
            ! First character? probably insert a separator.
            if (idestpos .gt. 0) then
              idestpos = idestpos + 1;
              sdest(idestpos:idestpos) = char(0)
            end if

            ! Transfer the characters
            idestpos = idestpos + 1
            sdest(idestpos:idestpos) = ssource(isourcepos:isourcepos)
            isourcepos = isourcepos+1

            idestpos = idestpos + 1
            sdest(idestpos:idestpos) = ssource(isourcepos:isourcepos)
            isourcepos = isourcepos+1
            
            ! We are now in two-character-expression mode
            ccurrentgroup = CHARGROUP_TWOCHAREXP
            
          end if
        end if
        
        if (j .eq. 0) exit
        
      end do
    
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine process_chargroup4 (ccurrentgroup,ssource,isourcepos,imaxlen,sdest,idestpos)
    
      ! Processes characters in character group 4.
      ! The destination and the pointers are modified according
      ! to the current state.
      
      ! Current character group
      integer, intent(inout) :: ccurrentgroup

      ! Source string
      character(len=*), intent(in) :: ssource
      
      ! Source pointer
      integer, intent(inout) :: isourcepos
      
      ! Length of ssource
      integer, intent(inout) :: imaxlen

      ! Destination string
      character(len=*), intent(inout) :: sdest
      
      ! Destination pointer
      integer, intent(inout) :: idestpos
      
      ! First character or new word (sequence)? probably insert a separator.
      if ((idestpos .gt. 1) .and. (ccurrentgroup .lt. CHARGROUP_INWORD)) then
        idestpos = idestpos + 1;
        sdest(idestpos:idestpos) = char(0)
      end if

      ! Is that an escape char?
      if (ssource(isourcepos:isourcepos) .eq. '\') then
        ! What about the next char?
        if (isourcepos .lt. imaxlen) then
          ! Take the next character as it is.
          idestpos = idestpos + 1
          sdest(idestpos:idestpos) = ssource(isourcepos+1:isourcepos+1)
          isourcepos = isourcepos+2
        else
          ! Ignore the character
          isourcepos = isourcepos+1
        end if
      else
        ! Transfer the character
        idestpos = idestpos + 1
        sdest(idestpos:idestpos) = ssource(isourcepos:isourcepos)
        
        isourcepos = isourcepos+1
      end if

      ! We are now in in-word mode.
      ccurrentgroup = CHARGROUP_INWORD
    
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine trigger_chargroup5 (ccurrentgroup,schars,btriggered,ssource,isourcepos,imaxlen,&
        sdest,idestpos,ccurrentquote,bdequote)
    
      ! Checks if the next token belongs to this character group.
      ! If that is the case, the destination and the pointers are modified according
      ! to the current state.
      
      ! Current character group
      integer, intent(inout) :: ccurrentgroup
      
      ! Characters in this group.
      character(len=*), intent(in) :: schars
      
      ! Returns if the group is triggered.
      logical, intent(out) :: btriggered
      
      ! Source string
      character(len=*), intent(in) :: ssource
      
      ! Source pointer
      integer, intent(inout) :: isourcepos
      
      ! Length of ssource
      integer, intent(inout) :: imaxlen

      ! Destination string
      character(len=*), intent(inout) :: sdest
      
      ! Destination pointer
      integer, intent(inout) :: idestpos
      
      ! Last found quotation character. =0 if no quotation active.
      character, intent(inout) :: ccurrentquote
      
      ! Determins whether to de-quote the string or not.
      logical, intent(in) :: bdequote
      
      ! local variables
      integer :: i
      
      ! Result is the current state of ccurrentquote
      btriggered = .false.
      
      ! Chech the current character for the next state.
      !
      ! Are we quoted?
      if (ccurrentgroup .ne. CHARGROUP_INQUOTE) then
        i = index(schars,ssource(isourcepos:isourcepos))

        if (i .ne. 0) then
          ! Quotation character. Start dequoting.
          ccurrentquote = ssource(isourcepos:isourcepos)
          
          ! First character or new word (sequence)? probably insert a separator
          ! for the characters that follow.
          if ((idestpos .gt. 1) .and. (ccurrentgroup .lt. CHARGROUP_INWORD)) then
            idestpos = idestpos + 1;
            sdest(idestpos:idestpos) = char(0)
          end if
          
          ! We are now in in-quote mode.
          ccurrentgroup = CHARGROUP_INQUOTE
          btriggered = .true.
          
          if (.not. bdequote) then
            ! Transfer the character
            idestpos = idestpos + 1
            sdest(idestpos:idestpos) = ssource(isourcepos:isourcepos)
          end if
          
          ! otherwise ignore the quote character.
          isourcepos = isourcepos+1
          
        end if
      else
        ! We are in 'quoted' mode. Any "\0" char was put to the string before.
        btriggered = .true.
        
        ! Is that an escape char?
        if (ssource(isourcepos:isourcepos) .eq. '\') then

          if (.not. bdequote) then
            ! Transfer this character as it is.
            idestpos = idestpos + 1
            sdest(idestpos:idestpos) = ssource(isourcepos:isourcepos)
          end if
        
          ! What about the next char?
          if (isourcepos .lt. imaxlen) then
            ! Take the next character as it is.
            idestpos = idestpos + 1
            sdest(idestpos:idestpos) = ssource(isourcepos+1:isourcepos+1)
            isourcepos = isourcepos+2
          else
            ! Ignore the character
            isourcepos = isourcepos+1
          end if
        else
          ! Stop quoting?
          if (ssource(isourcepos:isourcepos) .eq. ccurrentquote) then
            ! Return to 'inword' mode. Next whitespace triggers next word.
            ccurrentgroup = CHARGROUP_INWORD
            ccurrentquote = char(0)

            if (.not. bdequote) then
              ! Transfer the character
              idestpos = idestpos + 1
              sdest(idestpos:idestpos) = ssource(isourcepos:isourcepos)
            end if
            
            ! otherwise ignore the quote character.
            isourcepos = isourcepos+1
            
          else
            ! Transfer the character
            idestpos = idestpos + 1
            sdest(idestpos:idestpos) = ssource(isourcepos:isourcepos)
            
            isourcepos = isourcepos+1
          end if
        end if
      end if
    
    end subroutine

    ! -------------------------------------------------------------------------
    subroutine trigger_chargroup6 (ccurrentgroup,schars,btriggered,ssource,isourcepos,imaxlen)
    
      ! Checks if the next token belongs to this character group.
      ! If that is the case, the destination and the pointers are modified according
      ! to the current state.
      
      ! Current character group
      integer, intent(inout) :: ccurrentgroup

      ! Characters in this group.
      character(len=*), intent(in) :: schars
      
      ! Returns if the group is triggered.
      logical, intent(out) :: btriggered
      
      ! Source string
      character(len=*), intent(in) :: ssource
      
      ! Source pointer
      integer, intent(inout) :: isourcepos
      
      ! Length of ssource
      integer, intent(inout) :: imaxlen
      
      ! local variables
      integer :: i
      
      i = index(schars,ssource(isourcepos:isourcepos))
      btriggered = i .ne. 0
      
      if (btriggered) then
        ! Skip everything
        isourcepos = imaxlen
      end if
      
    end subroutine

  end subroutine
  
  ! ***************************************************************************

  subroutine cmdprs_nexttoken (ssource,istart,iend,ilength,bback)
  
  !<description>
    ! Determins the bounds of the next token.
  !</description>
  
  !<input>
    ! A source string. Must have been split with cmdprs_complexSplit.
    character(len=*), intent(in) :: ssource

    ! OPTIONAL: Trimmed length of the string. Speeds up the routine
    ! if specified.
    integer, intent(in), optional :: ilength
  !</input>
  
  !<inputoutput>
    ! On input: Start index of the previous token or =0, if the first
    ! token should be determined.
    ! On output: Start index of the next token or =0 if there is no more
    ! token.
    integer, intent(inout) :: istart
    
    ! On input: End index of the previous token or =0, if the first
    ! token should be determined.
    ! On output: End index of the next token or =0 if there is no more
    ! token.
    integer, intent(inout) :: iend
    
    ! OPTIONAL: Backward mode. Start from the end of the string.
    logical, intent(in), optional :: bback
  !</inputoutput>
  
  !</subroutine>
  
    integer :: ilen
    logical :: back
    if (present(ilength)) then
      ilen = ilength
    else
      ilen = len_trim(ssource)
    end if
    
    back = .false.
    if (present(bback)) back = bback
    
    if (.not. back) then
      
      ! Forward mode
      if (istart .eq. 0) then
        
        ! Find first token
        istart = 1
        iend = index (ssource,char(0))
        if (iend .eq. 0) then
          iend = istart
        else
          ! Go one back -- the \0.
          iend = iend - 1
        end if
        
      else
      
        ! Find next token
        istart = iend + 2
        if (istart .gt. ilen) then
          ! Finish
          istart = 0
          iend = 0
        else
          ! Find next end.
          iend = index (ssource(istart:),char(0))
          if (iend .eq. 0) then
            iend = ilength
          else
            ! Go one back -- the \0.
            iend = istart + iend - 2
          end if
        end if
      end if
    
    else

      ! Backward mode. Slightly easier due to well designed indices...
      
      if (iend .eq. 0) then
        
        ! Find first token
        iend = ilen
        istart = index (ssource,char(0),.true.) + 1
        
      else
      
        ! Find next token
        iend = istart - 2
        if (iend .lt. 1) then
          ! Finish
          istart = 0
          iend = 0
        else
          ! Find next end.
          istart = index (ssource(1:iend),char(0),.true.) + 1
        end if
      end if
    
    end if

  end subroutine

  ! ***************************************************************************

  !</subroutine>

  subroutine cmdprs_counttokens (ssource,ntokens)
  
  !<description>
    ! Calculates the number of tokens in ssource.
  !</description>
  
  !<input>
    ! A source string. Must have been split with cmdprs_complexSplit.
    character(len=*), intent(in) :: ssource
  !</input>
  
  !<output>
    ! Number of tokens. =0 if ssource does not contain text.
    integer, intent(out) :: ntokens
  !</output>
  
  !</subroutine>
  
    integer :: i
    
    ntokens = 0
    
    ! Return if there is no text.
    if (ssource .eq. "") return
  
    ! Count the number of \0 characters
    ntokens = 0
    do i = 1,len_trim(ssource)
      if (ssource(i:i) .eq. char(0)) then
        ntokens = ntokens + 1
      end if
    end do
    
    ! Add 1 for the last token.
    ntokens = ntokens + 1
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_commandmatch (ssource,ScommandTokens,bmatch)
  
  !<description>
    ! Checks if ssource matches a command.
    ! ssource must have been prepared with cmdprs_complexSplit.
    ! ScommandTokens is an array of strings which is compared to the
    ! substrings in ssource. A string "?" in ScommandTokens match any
    ! character sequence. A string "*" matches any number of arbitrary
    ! strings; remaining array entries are checked against the last tokens.
  !</description>
  
  !<input>
    ! A source string. Must have been split with cmdprs_complexSplit.
    character(len=*), intent(in) :: ssource
    
    ! Command subtokens to match against ssource.
    ! Examples:
    !   (/ "if", "(", "*", ")" /) matches an IF command.
    !   (/ "for", "(", "*" ")" /) matches a for command.
    ! Parsing stops with the first empty token in ScommandTokens.
    character(len=*), dimension(:), intent(in) :: ScommandTokens
  !</input>

  !<output>
    ! Whether the string matches or not.
    logical, intent(out) :: bmatch
  !</output>
  
    integer :: icurrenttokenstart, icurrenttokenend
    integer :: isubtoken1, isubtoken2, ilength
    
    bmatch = .true.
    isubtoken1 = 1
    
    ! Start to tokenise.
    icurrenttokenstart = 0
    icurrenttokenend = 0
    ilength = len_trim (ssource)
    
    do

      ! Get the next token    
      call cmdprs_nexttoken (ssource,icurrenttokenstart,icurrenttokenend,ilength)
      
      ! Check the token.
      !
      ! Was this the last token?
      if (ScommandTokens(isubtoken1) .eq. "") then
        if (icurrenttokenstart .eq. 0) then
          ! Done. String matches.
          return
        else
          ! There are more substrings than tokens. String does not match.
          bmatch = .false.
          return
        end if
      else if (icurrenttokenstart .eq. 0) then
        ! There are more tokens than substrings. String does not match.
        bmatch = .false.
        return
      end if
      
      ! Is this a placeholder for one word?
      if (ScommandTokens(isubtoken1) .eq. "?") then
        ! Go to next token.
        isubtoken1 = isubtoken1 + 1
        cycle
      end if
      
      ! Is this a placeholder for arbitrary many words?
      if (ScommandTokens(isubtoken1) .eq. "*") then
        ! Scan from the end.
        exit
      end if

      ! Does the token match?
      if (trim(ScommandTokens(isubtoken1)) .ne. &
          sys_upcase (ssource(icurrenttokenstart:icurrenttokenend))) then
        ! No! Cancel matching.
        bmatch = .false.
        return
      end if
      
      ! Next one.
      isubtoken1 = isubtoken1 + 1
      
    end do

    ! Find the last token to check
    do isubtoken2 = isubtoken1+1, ubound(ScommandTokens,1)
      if (ScommandTokens(isubtoken2) .eq. "") exit
    end do
    isubtoken2 = isubtoken2 - 1

    ! Start to tokenise from the end
    icurrenttokenend = 0
    icurrenttokenstart = 0
    
    do

      ! Get the next token    
      call cmdprs_nexttoken (ssource,icurrenttokenstart,icurrenttokenend,ilength,.true.)
      
      ! Is this a placeholder for one word?
      if (ScommandTokens(isubtoken2) .eq. "?") then
        ! Go to next token.
        isubtoken1 = isubtoken1 - 1
        cycle
      end if
      
      ! Is this a placeholder for arbitrary many words?
      if (ScommandTokens(isubtoken2) .eq. "*") then
        ! We reached the "*". String matches.
        return
      end if

      ! Does the token match?
      if (trim(ScommandTokens(isubtoken2)) .ne. &
          sys_upcase (ssource(icurrenttokenstart:icurrenttokenend))) then
        ! No! Cancel matching.
        bmatch = .false.
        return
      end if
      
      ! Next one.
      isubtoken2 = isubtoken2 - 1
      
    end do
  
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
      call cmdprs_complexSplit (stemp,sline,1,.true.)
      
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
    
    ! Command string. Splitted and dequoted.
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

    integer :: istart, iend, istart2, iend2, istart3, iend3, icmdIndex
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
    logical :: bmatch
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
    
    ! Name of the section corresponding to the current nesting level
    if (inestlevel .eq. 0) then
      ssectionname = ""
    else
      ssectionname = trim(sys_siL(inestlevel,10))
    end if

    ! Try to match the command
    call cmdprs_commandmatch (scommand,spatternINT,bmatch)
    
    if (bmatch) then
      ! Create INT variable.
      istart = 0
      iend = 0
      call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
      call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
      call collct_setvalue_int (rcmdStatus%rcollection, scommand(istart:iend), &
          0, .true., ssectionname=ssectionname) 
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommand,spatternDOUBLE,bmatch)
    
      if (bmatch) then
        ! Create DOUBLE variable.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call collct_setvalue_real (rcmdStatus%rcollection, scommand(istart:iend), &
            0.0_DP, .true., ssectionname=ssectionname) 
      end if
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommand,spatternSTRING,bmatch)
    
      if (bmatch) then
        ! Create STRING variable.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call collct_setvalue_string (rcmdStatus%rcollection, scommand(istart:iend), &
            "", .true., ssectionname=ssectionname) 
      end if
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommand,spatternINT2,bmatch)
    
      if (bmatch) then
        ! Create INT variable and assign.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call collct_setvalue_int (rcmdStatus%rcollection, scommand(istart:iend), &
            0, .true., ssectionname=ssectionname)
        istart2 = 0
        iend2 = 0
        call tpsym_evalExpression (scommand(istart:),rcmdStatus%rcollection,inestlevel,istart2,iend2,rvalue)
      end if
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommand,spatternDOUBLE2,bmatch)
    
      if (bmatch) then
        ! Create DOUBLE variable and assign.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call collct_setvalue_real (rcmdStatus%rcollection, scommand(istart:iend), &
            0.0_DP, .true., ssectionname=ssectionname) 
        istart2 = 0
        iend2 = 0
        call tpsym_evalExpression (scommand(istart:),rcmdStatus%rcollection,inestlevel,istart2,iend2,rvalue)
      end if
    end if

    if (.not. bmatch) then
      call cmdprs_commandmatch (scommand,spatternSTRING2,bmatch)
    
      if (bmatch) then
        ! Create STRING variable and assign.
        istart = 0
        iend = 0
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call cmdprs_nexttoken (scommand,istart,iend,len(scommand))
        call collct_setvalue_string (rcmdStatus%rcollection, scommand(istart:iend), &
            "", .true., ssectionname=ssectionname) 
        istart2 = 0
        iend2 = 0
        call tpsym_evalExpression (scommand(istart:),rcmdStatus%rcollection,inestlevel,istart2,iend2,rvalue)
      end if
    end if

    if (bworkflowAllowed .and. (rcmdStatus%ierror .eq. 0)) then
      
      ! 2nd chance...

      if (.not. bmatch) then
        
        call cmdprs_commandmatch (scommand,spatternBEGIN,bmatch)
      
        if (bmatch) then
          ! Get the subblock and execute
          call cmdprs_getNextBlock (rcmdblock,iline,istart2,iend2)
          
          ! Create new nest level in the collection
          call collct_addsection (rcmdStatus%rcollection, trim(sys_siL(inestlevel+1,10)))
          
          ! Execute
          call cmdprs_parsecmdblock (rcmdStatus,rcmdblock,inestlevel+1,istart2,iend2,rvalue)

          ! Remove local symbols
          call collct_deletesection (rcmdStatus%rcollection, trim(sys_siL(inestlevel+1,10)))

          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if
          
      if (.not. bmatch) then

        call cmdprs_commandmatch (scommand,spatternFOR,bmatch)
      
        if (bmatch) then
          ! For loop. Fetch the command block from the next lines.
          call cmdprs_getNextBlock (rcmdblock,iline+1,istart2,iend2)
          
          ! Process the command
          call cmdprs_doFor (rcmdStatus,rcmdblock,inestlevel,scommand,iline,istart2,iend2)
          
          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if
          
      if (.not. bmatch) then

        call cmdprs_commandmatch (scommand,spatternWHILE,bmatch)

        if (bmatch) then
          ! While loop. Fetch the command block from the next lines.
          call cmdprs_getNextBlock (rcmdblock,iline+1,istart2,iend2)

          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if

      if (.not. bmatch) then

        call cmdprs_commandmatch (scommand,spatternDO,bmatch)
      
        if (bmatch) then
          ! Do-while loop. Fetch the command block from the next lines.
          call cmdprs_getNextBlock (rcmdblock,iline+1,istart2,iend2)

          ! Move the current to the end of the block.
          iline = iend2
        end if
      end if

      if (.not. bmatch) then

        call cmdprs_commandmatch (scommand,spatternIF,bmatch)
      
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
      call cmdprs_splitlinedirect (scommand,p_Sargs)
      if (associated(p_Sargs)) then
        call cmdprs_dospecialcommand (rcmdStatus,p_Sargs,ierror)
        call cmdprs_releaseargs (p_Sargs)
        bmatch = ierror .eq. 0
        if (ierror .gt. 1) then
          ! =0: ok, =1: not found. >2: real error, pass to caller
          rcmdStatus%ierror = ierror
          bmatch = .true.
        end if
      end if
      
    end if
    
    if (.not. bmatch) then
      ! Try it as an immediate expression.
      ! Processes e.g. variable assignments etc...
      istart = 0
      iend = 0
      call tpsym_evalExpression (scommand,rcmdStatus%rcollection,inestlevel,istart,iend,rvalue)
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
      rcmdblock,inestlevel,iblockstart,iblockend,rreturnValue)
  
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
    
    ! Start line in the command block. If not specified, the first line is assumed.
    integer, intent(in), optional :: iblockstart
    
    ! End line in the command block. If not specified, the last line is assumed.
    integer, intent(in), optional :: iblockend
    
    ! Level of nesting
    integer, intent(in), optional :: inestlevel

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
        call cmdprs_parsecmdblock (rcmdStatus,rcmdblock,inestlevel,rreturnvalue=rreturnvalue)
  
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
    integer :: ierror

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
      
      if (len(sinput) .ne. 0) then
      
        ! Ignore comments
        if (getchar (sinput,1) .ne. "#") then

          ! Parse the line.
          call cmdprs_splitline (sinput,p_Sargs)
          
          ! A command involved?
          if (associated(p_Sargs)) then
          
            ! Execute the command
            call cmdprs_dospecialcommand (rcmdStatus,p_Sargs,ierror)
          
            ! Release memory
            call cmdprs_releaseargs(p_Sargs)
          end if
          
        end if
      
      end if
      
      call output_lbrk()
     
    end do
    
    ! Release memory.
    sinput = ""
  
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_parseIndexString1D (sstring,svecname,p_IvecComponents)
  
  !<description>
    ! Parses a string as 1D array. sstring may have the following form:
    ! "SPECIFIER" - returns svecname="SPECIFIER" and p_IvecComponents => null().
    ! "SPECIFIER[i]" - returns svecname="SPECIFIER" and p_IvecComponents => [i,1].
    ! "SPECIFIER[i,j,...]" - returns svecname="SPECIFIER" and p_IvecComponents => [i,j,...].
    ! "SPECIFIER[i:j]" - returns svecname="SPECIFIER" and p_IvecComponents => [i,...,j].
  !</description>

  !<input>
    ! Input string.
    character(len=SYS_STRLEN), intent(in) :: sstring
  !</input>

  !<output>
    ! The specifier in the string. ="" if sstring is invalid.
    character(len=SYS_STRLEN), intent(out) :: svecname
    
    ! Pointer to an array collecting the indices or NULL if no indices are specified.
    ! Must be released by the caller!
    integer, dimension(:), pointer :: p_IvecComponents
  !</output>

    integer :: i,j,k,i1,i2,icount
    character(len=SYS_STRLEN) :: sstr,ssubstr

    ! Reset the output.
    svecname = ""
    nullify(p_IvecComponents)
    
    ! Find a bracket.
    i = index(sstring,"[")
    if (i .eq. 0) then
      ! Only the name is given.
      svecname = trim(adjustl(sstring))
      return
    else
      ! Get the name and the specifier. Remove all spaces from the specifier.
      svecname = trim(adjustl(sstring(1:i-1)))
      sstr = ""
      k = 0
      do j = i,len_trim(sstring)
        select case (ichar(sstring(j:j)))
        case (ichar(' '),ichar('['),ichar(']'))
          ! Do nothing.
        case (ichar('A'):ichar('Z'),ichar('_'),ichar('0'):ichar('9'),ichar(':'))
          ! Copy
          k=k+1
          sstr(k:k) = sstring(j:j)
        case (ichar(','))
          ! Replace by space
          k = k+1
          sstr(k:k) = ' '
        case default
          ! Stop here, the string is invalid!
          svecname = ""
          return
        end select
      end do
      
      ! Ok. Now find the number sets. Orient at the whitespaces.
      i = 1
      icount = 0
      do while (i .le. k)
        
        ! Get the first substring.
        read(sstr,'(A)') ssubstr
        
        ! ":" included? Replace by ' ' and read in the set.
        j = index(ssubstr,":")
        if (j .ne. 0) then
          ssubstr(j:j) = "  "
          read(ssubstr,*) i1,i2
          icount = icount + abs(i2 - i1) + 1
        else
          ! Only one entry
          icount = icount + 1
        end if
        
        ! Next token
        i = i + len_trim(ssubstr)
        
      end do
      
      ! Allocate memory for the numbers
      allocate (p_IvecComponents(icount))
      
      ! Get the numbers.
      icount = 0
      i = 1
      do while (i .le. k)
        
        ! Get the first substring.
        read(sstr,'(A)') ssubstr
        
        ! ":" included? Replace by ' ' and read in the set.
        j = index(ssubstr,":")
        if (j .ne. 0) then
        
          ssubstr(j:j) = " "
          read(ssubstr,*) i1,i2
          
          if (i1 .le. i2) then
            do j = i1,i2
              icount = icount + 1
              p_IvecComponents(icount) = j
            end do
          else
            do j = i2,i1,-1
              icount = icount + 1
              p_IvecComponents(icount) = j
            end do
          end if
        else
          ! Only one entry
          icount = icount + 1
          read (ssubstr,*) p_IvecComponents(icount)
        end if
        
        ! Next token
        i = i + len_trim(ssubstr)
        
      end do
      
    end if
    
  end subroutine

  ! ***************************************************************************

  recursive subroutine cmdprs_dospecialcommand (rcmdStatus,Sargs,ierror)
  
  !<description>
    ! Executes a command.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</inputoutput>
  
  !<output>
    ! Error flag.
    ! =0: no error.
    ! =1: Command not found
    integer, intent(out) :: ierror
  !</output>

    type(VARYING_STRING) :: scmd, sarg
    
    ierror = 0
    
    ! Get the command
    scmd = Sargs(1)
    call TOUPPER(scmd)
    
    if (scmd .eq. "EXIT") then
      rcmdStatus%bterminate = .true.
      return
    end if
    
    if (scmd .eq. "PRINT") then
      if (size (Sargs) .gt. 1) then
        sarg = cmdprs_dequote(Sargs(2))
        call output_line (char(sarg))
        sarg = ""
      else
        call output_lbrk ()
      end if
      return
    end if

    if (scmd .eq. "MEMINFO") then
      call output_lbrk ()
      call storage_info(.true.)
      return
    end if

    if (scmd .eq. "SET") then
      call cmdprs_do_set (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "ECHO") then
      if (size (Sargs) .gt. 1) then
        if (sys_upcase(char(Sargs(2))) .eq. "ON") then
          rcmdStatus%becho = .true.
        else if (sys_upcase(char(Sargs(2))) .eq. "OFF") then
          rcmdStatus%becho = .false.
        end if
      else
        if (rcmdStatus%becho) then
          call output_line ("Echo is on.")
        else
          call output_line ("Echo is off.")
        end if
      end if
      return
    end if

    if (scmd .eq. "RUN") then
      call cmdprs_do_run (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "SHOW") then
      call cmdprs_do_show (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "INFO") then
      call cmdprs_do_info (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "DELETE") then
      call cmdprs_do_delete (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "DESTROY") then
      call cmdprs_do_destroy (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "HELP") then
      call cmdprs_do_help (Sargs)
      return
    end if

    if (scmd .eq. "READ2DPRM") then
      call cmdprs_do_read2dprm (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "READ2DTRI") then
      call cmdprs_do_read2dtri (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "MESHREFINE") then
      call cmdprs_do_meshrefine (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "MESHHIERARCHY") then
      call cmdprs_do_meshhierarchy (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "FESPACE") then
      call cmdprs_do_fespace (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "FEHIERARCHY") then
      call cmdprs_do_fehierarchy (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "MLEVELPRJHIERARCHY") then
      call cmdprs_do_mlevelprjhierarchy (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "READBLOCKVECTOR") then
      call cmdprs_do_readblockvector (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "WRITEBLOCKVECTOR") then
      call cmdprs_do_writeblockvector (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "CREATEBLOCKVECTOR") then
      call cmdprs_do_createblockvector (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "COPYVECTOR") then
      call cmdprs_do_copyvector (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "INTERPOLATEVECTOR") then
      call cmdprs_do_interpolatevector (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "L2PROJECTION") then
      call cmdprs_do_l2projection (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "WRITEUCD") then
      call cmdprs_do_writeucd (rcmdStatus,Sargs)
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
  
  recursive subroutine tpsym_evalExpression (sstring,rcollection,inestlevel,istart,iend,rvalue)
  
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
    ! Collection containing symbols.
    type(t_collection), intent(inout) :: rcollection

    ! Level of nesting
    integer, intent(in) :: inestlevel
  !</input>
  
  !<input
  
  !<output>
    ! The value of the expression.
    type(t_symbolValue), intent(out) :: rvalue
  !</output>
  
  !</subroutine>
  
    ! local variables
    integer :: iprevopstart,iprevopend,icurropstart,icurropend,isymbolstart,isymbolend
    type(t_symbolValue) :: rvaluetemp
    integer :: ileftprec, irightprec, ilength
    
    ilength = len_trim(sstring)
    
    ! Start with the first token?
    if (istart .eq. 0) then
      call cmdprs_nexttoken (sstring,istart,iend,ilength)
    end if

    ! Current token.
    isymbolstart = istart
    isymbolend = iend
    
    ! Get the previous operator.
    iprevopstart = istart
    iprevopend = iend
    call cmdprs_nexttoken (sstring,iprevopstart,iprevopend,ilength,.true.)
    
    if (sstring(istart:iend) .eq. "(") then
      ! Evaluate the value in the brackets.
      call tpsym_evalExpressionBrackets (sstring,rcollection,inestlevel,istart,iend,rvalue)
      ! istart/iend now points to the last token -- the closing bracket.
    else
      ! Parse the symbol, create a value.
      call tpsym_parseSymbol (sstring(isymbolstart:isymbolend),rcollection,inestlevel,rvalue)
    end if
    if (rvalue%ctype .eq. STYPE_INVALID) return
    
    ! Repeat to evaluate all following operators.
    do while (istart .ne. 0)
      ! We have a value. Get the next operator.
      call cmdprs_nexttoken (sstring,istart,iend,ilength)
      icurropstart = istart
      icurropend = iend
      
      ! What is more important, the left or the rigt operator?
      !
      ! Left more important -> Stop here, return, value will be used.
      ! Right more important -> Recursively evaluate the expression on the right
      !    and evaluate the value of the operator.
      if (iprevopstart .gt. 0) then
        if (icurropstart .eq. 0) return
        if (getOperatorWeight(sstring(iprevopstart:iprevopend)) .gt. &
            getOperatorWeight(sstring(icurropstart:icurropend))) return
      end if
      
      ! Check for immediate operators.
      if (sstring(icurropstart:icurropend) .eq. "++") then
      
        ! Increase variable, return old value.
        rvaluetemp = 1
        rvalue = rvalue + rvaluetemp
        call tpsym_saveSymbol (rvalue,inestlevel,rcollection)
        rvalue = rvalue - rvaluetemp
        
        ! Next token
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
        if (istart .eq. 0) return
        
      elseif (sstring(icurropstart:icurropend) .eq. "--") then
      
        ! Increase variable, return old value.
        rvaluetemp = 1
        rvalue = rvalue - rvaluetemp
        call tpsym_saveSymbol (rvalue,inestlevel,rcollection)
        rvalue = rvalue + rvaluetemp

        ! Next token
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
        if (istart .eq. 0) return
        
      else      

        ! Right operator more important. Evaluate the right expression first and
        ! calculate a new value with the operator.
        call cmdprs_nexttoken (sstring,istart,iend,ilength)
        if (istart .eq. 0) return
        
        call tpsym_evalExpression (sstring,rcollection,inestlevel,istart,iend,rvaluetemp)
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
          call tpsym_saveSymbol (rvalue,inestlevel,rcollection)
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
      else if ((soperator .eq. "+") .or. (soperator .eq. "-") .or. &
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
  
  recursive subroutine tpsym_evalExpressionBrackets (sstring,rcollection,inestlevel,istart,iend,rvalue)
  
  !<description>
    ! Parses an expression in brackets.
  !</description>

  !<input>
    ! String containing symbols. The string must have been split with cmdprs_complexSplit
    ! and must completely be trimmed!
    character(len=*), intent(in) :: sstring
  
    ! Collection containing symbols.
    type(t_collection), intent(inout) :: rcollection

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
    integer :: istart2,iend2
    integer :: ibracket,ilength
    
    ! Get basic information
    ilength = len(sstring)
    
    ! Expression in brackets. Most important.
    !
    ! Search closing bracket -- or inner expression(s) to be evaluated.
    istart2 = istart
    iend2 = iend
    ibracket = 1
    call cmdprs_nexttoken (sstring,istart2,iend2,ilength,.true.)
    do while ((ibracket .gt. 0) .and. (istart2 .ne. 0))
      if (sstring(istart2:istart2) .eq. "(") then
        ibracket = ibracket + 1
      else if (sstring(istart2:istart2) .eq. ")") then
        ibracket = ibracket - 1
        if (ibracket .eq. 0) exit
      end if
      
      ! Ignore other tokens.
      call cmdprs_nexttoken (sstring,istart2,iend2,ilength,.true.)
    end do
    
    ! Return if ibracket <> 0, expression invalid.
    if (ibracket .gt. 0) return
    
    ! Evaluate recursively.
    istart = 0
    iend = 0
    call tpsym_evalExpression (sstring(iend+2:istart2-2),rcollection,inestlevel,&
        istart,iend,rvalue)
    
    ! Return the position of the closing bracket.
    istart = istart2
    iend = iend2

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

    ! Command string. Splitted and dequoted.
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
          call cmdprs_parsecmdblock (rcmdStatus,rcmdblock,inestlevel,&
              iblockstart,iblockend,rreturnvalue)

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

  ! ***************************************************************************

  !<subroutine>

  subroutine cmdprs_do_help (Sargs)
  
  !<description>
    ! Command: HELP.
  !</description>
  
  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
  !</subroutine>  
  
    character(len=20) :: sargformatted
  
    call output_lbrk ()

    if (size(Sargs) .eq. 1) then
      call output_line ("FEAT2 command parser.")
      call output_line ("---------------------")
      call output_line ("The following commands are available. For a specific help")
      call output_line ("type ""HELP [command]"".")
      call output_lbrk ()
      call output_line ("  help               - This help page.")
      call output_line ("  exit               - Exits the command line.")
      call output_line ("  meminfo            - Information about memory management.")
      call output_line ("  run                - Execute script.")
      call output_line ("  print              - Print some text.")
      call output_line ("  set                - Modify/set/print environment.")
      call output_line ("  show               - Show environment variable (if possible).")
      call output_line ("  delete             - Delete variable from environment.")
      call output_line ("  destroy            - Destroys an environment structure (extended delete).")
      call output_line ("  read2dprm          - Read 2D .prm file.")
      call output_line ("  read2dtri          - Read 2D .tri file.")
      call output_line ("  meshrefine         - Refine a mesh.")
      call output_line ("  meshhierarchy      - Create a mesh hierarchy.")
      call output_line ("  fespace            - Create a FE space.")
      call output_line ("  fehierarchy        - Create a FE hierarchy.")
      call output_line ("  mlevelprjhierarchy - Create a multilevel projection hierarchy.")
      call output_line ("  createblockvector  - Create an empty block vector.")
      call output_line ("  readblockvector    - Read a block vector from a file.")
      call output_line ("  writeblockvector   - Write a block vector to a file.")
      call output_line ("  copyvector         - Copy a vector to another.")
      call output_line ("  interpolatevector  - Interpolate a vector to another level.")
      call output_line ("  l2projection       - Appies an L2 projection to a vector.")
      call output_line ("  writeucd           - Writes a postprocessing file.")
    
    else
    
      call cmdprs_getparam (Sargs,2,sargformatted,.true.,.false.)
      
      if (sargformatted .eq. "EXIT") then
        call output_line ("EXIT - Close interactive command line, return to terminal.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  EXIT")
      
        return
      end if

      if (sargformatted .eq. "RUN") then
        call output_line ("RUN - Execute a script on hard disc.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  RUN ""[filename/path]""")
      
        return
      end if

      if (sargformatted .eq. "PRINT") then
        call output_line ("PRINT - prints some text to the terminal.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  PRINT ""[some text]""")
      
        return
      end if

      if (sargformatted .eq. "SET") then
        call output_line ("SET - Modify or show environment.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("* To show statistics about the environment:")
        call output_line ("      SET")
        call output_lbrk ()
        call output_line ("*To set/modify an integer variable:")
        call output_line ("      SET INT [name] [value]")
        call output_lbrk ()
        call output_line ("*To set/modify an double precision variable:")
        call output_line ("      SET DOUBLE [name] [value]")
        call output_lbrk ()
        call output_line ("*To set/modify a string variable:")
        call output_line ("      SET STRING [name] [value]")
      
        return
      end if
    
      if (sargformatted .eq. "DELETE") then
        call output_line ("DELETE - Delete a variable from the environment.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  DELETE [name]")
        call output_lbrk ()
        call output_line ("WARNING: This does not release probably associated memory!")
      
        return
      end if
    
      if (sargformatted .eq. "DESTROY") then
        call output_line ("DESTROY - Releases memory associated to an environment variable.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  DESTROY [name]")
        call output_lbrk ()
        call output_line ("Automatically detects the type of an environment variable and")
        call output_line ("releases the associated structure from memory (if there is one).")
        call output_line ("The variable is removed from the environment")
      
        return
      end if
    
      if (sargformatted .eq. "SHOW") then
        call output_line ("SHOW - Show content about an environment vairbale.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  SHOW [name]")
        call output_lbrk ()
        call output_line ("If possible, this routine prints information about an")
        call output_line ("environment variable to the terminal. For standard types")
        call output_line ("(int, double, string,...), the value is returned.")
        call output_line ("For extended types (e.g. meshes), status information")
        call output_line ("is returned.")
      
        return
      end if
    
      if (sargformatted .eq. "INFO") then
        call output_line ("INFO - Detailed information about environment vairbales.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  SHOW [name]")
        call output_lbrk ()
        call output_line ("Shows detailed information about an environment variable.")
        call output_line ("This includes type, value,...")
      
        return
      end if
    
      if (sargformatted .eq. "READ2DPRM") then
        call output_line ("READ2DPRM - Read 2D .PRM file.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  READ2DPRM [varname] [filename]")
        call output_lbrk ()
        call output_line ("Read in a .prm file and store it using the name [variable].")
      
        return
      end if

      if (sargformatted .eq. "READ2DTRI") then
        call output_line ("READ2DTRI - Read 2D .TRI file.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("  READ2DTRI [varname] [filename] [ BOUNDARY varbd ]")
        call output_lbrk ()
        call output_line ("Read in a .tri file and store it using the name [variable].")
        call output_line ("If BOUNDARY is specified, the triangulation is connected")
        call output_line ("to a boundary object varbd.")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("  read2dprm myprm ""myprmfile.prm""")
        call output_line ("  read2dtri mytri ""mytrifile.tri"" boundary myprm")
      
        return
      end if

      if (sargformatted .eq. "MESHREFINE") then
        call output_line ("MESHREFINE - Refine a mesh.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   meshrefine [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Refine a mesh with a given method one or multiple times.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --TIMES [varname] ...")
        call output_line ("      Refine the mesh [varname] times. Default is ""--TIMES 1"".")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use the specified boundary object [varname] fdor boundary refinement.")
        call output_line ("      If not specified, no boundary object is used.")
        call output_lbrk ()
        call output_line ("  ... --METHOD [varname] ...")
        call output_line ("      Use a specific method. Default is ""--METHOD 2LEVELORDERED"". Possible choices:")
        call output_line ("      2LEVELORDERED   - Use 2-level-ordering refinement.")
      
        return
      end if

      if (sargformatted .eq. "MESHHIERARCHY") then
        call output_line ("MESHHIERARCHY - Create a mesh hierarchy.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   meshhierarchy [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Refine a mesh with a given method one or multiple times.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --MESH [varname] ...")
        call output_line ("      Use mesh [varname] as coarse mesh of teh hierarchy.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use boundary definition [varname] for refinements.")
        call output_lbrk ()
        call output_line ("  ... --LEVELS [varname] ...")
        call output_line ("      Creates a hierarchy of [varname] levels. Default is ""--LEVELS 1"".")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use the specified boundary object [varname] fdor boundary refinement.")
        call output_line ("      If not specified, no boundary object is used.")
        call output_lbrk ()
        call output_line ("  ... --METHOD [varname] ...")
        call output_line ("      Use a specific method. Default is ""--METHOD 2LEVELORDERED"". Possible choices:")
        call output_line ("      2LEVELORDERED   - Use 2-level-ordering refinement.")
      
        return
      end if

      if (sargformatted .eq. "FESPACE") then
        call output_line ("FESPACE - Create an FE space.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   fespace [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Creates an FE space that can be used to set up matrices/vectors.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --MESH [varname] ...")
        call output_line ("      Use mesh [varname] as coarse mesh of the hierarchy.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use boundary definition [varname] for refinements.")
        call output_lbrk ()
        call output_line ("  ... --COMPONENTS [varname] ...")
        call output_line ("      Creates a FE space with [varname] components. Mandatory argument.")
        call output_line ("      Must be defined in advance to all element/cubature specifications.")
        call output_lbrk ()
        call output_line ("  ... --TRIELEMENTS EL_xxx EL_xxx EL_xxx ...")
        call output_line ("      Specifies a list of TRI/TETRA elemnent types to be used for")
        call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --QUADELEMENTS EL_xxx EL_xxx EL_xxx ...")
        call output_line ("      Specifies a list of QUAD/HEXA elemnent types to be used for")
        call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --TRICUB CUB_xxx CUB_xxx CUB_xxx ...")
        call output_line ("      Specifies a list of TRI/TETRA cubature formulas to be used for")
        call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
        call output_line ("      If not defined, the default cubature formula is used.")
        call output_lbrk ()
        call output_line ("  ... --QUADCUB CUB_xxx CUB_xxx CUB_xxx ...")
        call output_line ("      Specifies a list of QUAD/HEXA cubature formulas to be used for")
        call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
        call output_line ("      If not defined, the default cubature formula is used.")
      
        return
      end if

      if (sargformatted .eq. "FEHIERARCHY") then
        call output_line ("FEHIERARCHY - Create an FE hierarchy.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   fehierarchy [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Creates an FE space that can be used to set up matrices/vectors.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --MESHHIERARCHY [varname] ...")
        call output_line ("      Use mesh hierarchy [varname] as coarse mesh of the hierarchy.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use boundary definition [varname] for refinements.")
        call output_lbrk ()
        call output_line ("  ... --COMPONENTS [varname] ...")
        call output_line ("      Creates a FE space with [varname] components. Mandatory argument.")
        call output_line ("      Must be defined in advance to all element/cubature specifications.")
        call output_lbrk ()
        call output_line ("  ... --TRIELEMENT EL_xxx EL_xxx EL_xxx ...")
        call output_line ("      Specifies a list of TRI/TETRA elemnent types to be used for")
        call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --QUADELEMENT EL_xxx EL_xxx EL_xxx ...")
        call output_line ("      Specifies a list of QUAD/HEXA elemnent types to be used for")
        call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --TRICUB CUB_xxx CUB_xxx CUB_xxx ...")
        call output_line ("      Specifies a list of TRI/TETRA cubature formulas to be used for")
        call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --QUADCUB CUB_xxx CUB_xxx CUB_xxx ...")
        call output_line ("      Specifies a list of QUAD/HEXA cubature formulas to be used for")
        call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
      
        return
      end if

      if (sargformatted .eq. "MLEVELPRJHIERARCHY") then
        call output_line ("MLEVELPRJHIERARCHY - Create an multilevel projection hierarchy.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   mlevelprjhierarchy [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Creates a multilevel projection hierarchy that can be used to transfer")
        call output_line ("soltion/rhs vectors from one level to another.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --FEHIERARCHY [varname] ...")
        call output_line ("      Create a multilevel projection hierarchy based on the")
        call output_line ("      FE space hierarchy [varname]. Mandatory argument")
      
        return
      end if

      if (sargformatted .eq. "CREATEBLOCKVECTOR") then
        call output_line ("CREATEBLOCKVECTOR - Create an empty block vector.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   createblockvector [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Creates an empty block vector..")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --FEHIERARCHY [varname] ...")
        call output_line ("      Defines the FE hierarchy, the vector is based on.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --LEVEL [varname] ...")
        call output_line ("      Level in the hierarchy specifying the discretisation of the vector.")
        call output_line ("      Default is max. level.")
      
        return
      end if

      if (sargformatted .eq. "READBLOCKVECTOR") then
        call output_line ("READBLOCKVECTOR - Read block vector from file.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   readblockvector [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Read a block vector from a file.")
        call output_line ("[varname] identifies the variable where to read data to.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --FILENAME [varname] ...")
        call output_line ("      Defines the filename of the file to read.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --UNFORMATTED ...")
        call output_line ("      Read a binary file. Default is formatted, human readable file.")
      
        return
      end if

      if (sargformatted .eq. "WRITEBLOCKVECTOR") then
        call output_line ("WRITEBLOCKVECTOR - Write block vector to file.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   writeblockvector [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Write a block vector to a file.")
        call output_line ("[varname] identifies the vector.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --FILENAME [varname] ...")
        call output_line ("      Defines the filename of the file to read.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --FORMAT [varname] ...")
        call output_line ("      Defines the format of the numbers in the file.")
        call output_line ("      Applies only if --UNFORMATTED is not specified.")
        call output_line ("      Default: E20.10")
        call output_lbrk ()
        call output_line ("  ... --UNFORMATTED ...")
        call output_line ("      Read a binary file. Default is formatted, human readable file.")
      
        return
      end if

      if (sargformatted .eq. "COPYVECTOR") then
        call output_line ("COPYVECTOR - Copy a vector.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   copyvector [varsource] [vardest]")
        call output_lbrk ()
        call output_line ("Copies vector [varsource] to vector [vardest].")
        call output_line ("The vector may be a full block vector or a single.")
        call output_line ("subvector.")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("    copyvector source dest")
        call output_line ("    copyvector source[2] dest[1]")
      
        return
      end if

      if (sargformatted .eq. "INTERPOLATEVECTOR") then
        call output_line ("INTERPOLATEVECTOR - Interpolate a vector.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   interpolatevector [varsource] [lvlsource] [vardest] [lvldest] [...options...]")
        call output_lbrk ()
        call output_line ("Interpolates vector [varsource] from level [lvlsource] to")
        call output_line ("level [lvldest] and writes the result to [vardest].")
        call output_line ("[mlprj] specifies the projection hierarchy to be used.")
        call output_line ("The method uses prolongation/interpolation for the level change.")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("    interpolatevector source 3 dest 5 --MLPRJHIERARCHY mhier --FEHIERARCHY feh")
        call output_lbrk ()
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --MLPRJHIERARCHY [varname] ...")
        call output_line ("      Defines the underlying projection hierarchy.")
        call output_line ("      Mandatory argument.")
        call output_lbrk ()
        call output_line ("  ... --FEHIERARCHY [varname] ...")
        call output_line ("      Defines the underlying FE hierarchy.")
        call output_line ("      Mandatory argument if source and destination level differ")
        call output_line ("      by more than one level.")
      
        return
      end if

      if (sargformatted .eq. "L2PROJECTION") then
        call output_line ("L2PROJECTION - Applies an L2 projection to a vector.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   l2projection [varsource] [vardest] [...options...]")
        call output_lbrk ()
        call output_line ("Interpolates vector [varsource] to [vardest] using an")
        call output_line ("L2-projection.")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("    interpolatevector source dest")
        call output_lbrk ()
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --VERBOSE ...")
        call output_line ("      Activate verbose output.")
        call output_lbrk ()
        call output_line ("  ... --RELERROR [error] ...")
        call output_line ("      Defines the relative accuracy of the projection.")
      
        return
      end if

      if (sargformatted .eq. "WRITEUCD") then
        call output_line ("WRITEUCD - Write a postprocessing file.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   writeucd [filename] [mesh] [type] [...options...]")
        call output_lbrk ()
        call output_line ("Writes a postproceessing file [filename] based on the")
        call output_line ("mesh [mesh]. [type] specifies the type of the output.")
        call output_line ("[mesh] may be either the name of a mesh or in an indexed")
        call output_line ("form a mesh in a mesh hierarchy, e.g. ""myhier[3]"".")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("    writeucd mypostproc.vtk mymesh vtk")
        call output_lbrk ()
        call output_line ("The following postprocessing output is possible:")
        call output_lbrk ()
        call output_line ("  [type] = ""vtk"": VTK-output")
        call output_line ("           ""gmv"": GMV output")
        call output_lbrk ()
        call output_line ("The following parameters are available in [...options...]:")
        call output_lbrk ()
        call output_line ("  --POINTDATASCALAR [name] [vector]")
        call output_line ("          Writes out subvector of vector [vector]")
        call output_line ("          as scalar array with the name [name]. May be specified")
        call output_line ("          more than once. The data is interpreted in the corner")
        call output_line ("          points of the elements. Example:")
        call output_line ("              --pointdatascalar vel_y myvec[2] 2")
        call output_lbrk ()
        call output_line ("  --POINTDATAVEC [name] [vector]")
        call output_line ("          Writes out a set of subvectors from vector [vector]")
        call output_line ("          as a vector field with name [name]. [subveccount] specifies")
        call output_line ("          the number of subvectors and [subvectors...] is a list")
        call output_line ("          of subvectors. Example to write out component 1 and 2:")
        call output_line ("              --pointdatavec velocity myvec[1,2]")
        call output_line ("          Example to write out component 1 till 3:")
        call output_line ("              --pointdatavec velocity myvec[1..3]")
        call output_lbrk ()
        call output_line ("  --CELLDATASCALAR [name] [vector] [subvector]")
        call output_line ("          Writes out subvector [subvector] of vector [vector]")
        call output_line ("          as scalar array with the name [name]. May be specified")
        call output_line ("          more than once. The data is interpreted in the elements.")
        call output_line ("          Example:")
        call output_line ("              --celldatascalar pressure myvec 3")
      
        return
      end if

      call output_line ("No specific help available.")
    end if
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine cmdprs_do_run (rcmdStatus,Sargs)
  
  !<description>
    ! Command: RUN.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    integer :: iunit
    logical :: bexists
    character(len=SYS_STRLEN) :: svalue
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the script name
    call cmdprs_getparam (Sargs,2,svalue,.false.,.true.)

    ! Open the file and read.
    inquire(file=trim(svalue), exist=bexists)
    if (.not. bexists) then
      call output_line ("File not found!")
    else
      
      ! Open
      call io_openFileForReading(svalue, iunit, .true.)

      ! Interpret the stream
      call cmdprs_parsestream (rcmdStatus,iunit,0)
      
      close (iunit)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_set (rcmdStatus,Sargs)
  
  !<description>
    ! Command: SET.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    integer :: iarg
    real(DP) :: darg
    character(len=COLLCT_MLNAME) :: svalname
    character(len=20) :: sformat
    character(len=SYS_STRLEN) :: svalue
  
    if (size(Sargs) .eq. 1) then
      call output_lbrk ()
      call output_line ("Environment statistics:")
      call output_lbrk ()
      ! Print the current status of the collection.
      call collct_printStatistics (rcmdStatus%rcollection)
      
    else
    
      if (size(Sargs) .lt. 4) then
        call output_line ("Not enough arguments.")
        return
      end if
      
      ! Get the next token(s)
      call cmdprs_getparam (Sargs,2,sformat,.true.,.false.)
      call cmdprs_getparam (Sargs,3,svalname,.false.,.true.)
      
      if (sformat .eq. "INT") then
        call cmdprs_getparam (Sargs,4,svalue,.false.,.true.)
        read (svalue,*) iarg
        call collct_setvalue_int (rcmdStatus%rcollection, &
            svalname, iarg, .true.)
        return
      end if
    
      if (sformat .eq. "DOUBLE") then
        call cmdprs_getparam (Sargs,4,svalue,.false.,.true.)
        read (svalue,*) darg
        call collct_setvalue_real (rcmdStatus%rcollection, &
            svalname, darg, .true.)
        return
      end if

      if (sformat .eq. "STRING") then
        call cmdprs_getparam (Sargs,4,svalue,.true.,.false.)
        call collct_setvalue_string (rcmdStatus%rcollection, &
            svalname, trim(svalue), .true.)
        return
      end if
    
    end if
  
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_delete (rcmdStatus,Sargs)
  
  !<description>
    ! Command: SET.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: svalname
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the next token(s)
    call cmdprs_getparam (Sargs,2,svalname,.true.,.false.)
    
    call collct_deletevalue (rcmdStatus%rcollection, svalname)
    
  end subroutine

  ! ***************************************************************************

  subroutine do_destroy (rcmdStatus,svalname,bverbose)
  
  !<description>
    ! Command handler: DESTROY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Name of the variable to destroy
    character(len=*), intent(in) :: svalname
    
    ! Print vebose messages
    logical, intent(in) :: bverbose
  !</input>
  
    ! local variables
    integer :: ctype
    logical :: bexists
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_interlevelProjectionHier), pointer :: p_rinterlevelProjectionHier
    type(t_vectorBlock), pointer :: p_rvectorBlock
  
    ! That's an awful task. Figure out the type at first.
    ! Then get the value and/or print information.
    ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists)
    
    if (.not. bexists) then
      if (bverbose) then
        call output_line ("Unknown environment variable!")
      end if
    else
      select case (ctype)
      case (COLLCT_INTEGER,COLLCT_REAL,COLLCT_STRING)
        ! Primitive type
        call collct_deletevalue (rcmdStatus%rcollection, svalname)
          
      case (COLLCT_BOUNDARY)
        ! Boundary object
        p_rboundary => collct_getvalue_bdry (rcmdStatus%rcollection, svalname)
        call boundary_release(p_rboundary)
        deallocate(p_rboundary)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)
          
      case (COLLCT_TRIA)
        ! Boundary object
        p_rtriangulation => collct_getvalue_tria (rcmdStatus%rcollection, svalname)
        call tria_done(p_rtriangulation)
        deallocate(p_rtriangulation)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_MSHHIERARCHY)
        ! Boundary object
        p_rmeshhierarchy => collct_getvalue_mshh (rcmdStatus%rcollection, svalname)
        call mshh_releaseHierarchy(p_rmeshhierarchy)
        deallocate(p_rmeshhierarchy)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_FESPACE)
        ! Boundary object
        p_rfeSpace => collct_getvalue_fesp (rcmdStatus%rcollection, svalname)
        call fesph_releaseFEspace(p_rfeSpace)
        deallocate(p_rfeSpace)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_FEHIERARCHY)
        ! Boundary object
        p_rfeHierarchy => collct_getvalue_feh (rcmdStatus%rcollection, svalname)
        call fesph_releaseHierarchy(p_rfeHierarchy)
        deallocate(p_rfeHierarchy)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_MLPRJHIERARCHY)
        ! Multilevel projection hierarchy
        p_rinterlevelProjectionHier => collct_getvalue_mlprjh (rcmdStatus%rcollection, svalname)
        call mlprj_releasePrjHierarchy(p_rinterlevelProjectionHier)
        deallocate(p_rinterlevelProjectionHier)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_BLKVECTOR)
        ! Multilevel projection hierarchy
        p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svalname)
        call lsysbl_releaseVector(p_rvectorBlock)
        deallocate(p_rvectorBlock)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case default
        call output_line ("Unknown type, variable cannot be destroyed!")         
      end select
      
    end if
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_destroy (rcmdStatus,Sargs)
  
  !<description>
    ! Command: DESTROY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: svalname
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the next token(s)
    call cmdprs_getparam (Sargs,2,svalname,.true.,.false.)
    
    call do_destroy (rcmdStatus,svalname,.true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_show (rcmdStatus,Sargs)
  
  !<description>
    ! Command: SHOW.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: svalname
    character(len=SYS_STRLEN) :: svalue
    integer :: ctype
    logical :: bexists
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the next token(s)
    call cmdprs_getparam (Sargs,2,svalname,.true.,.false.)
    
    ! That's an awful task. Figure out the type at first.
    ! Then get the value and/or print information.
    ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists)
    
    if (.not. bexists) then
      call output_line ("Unknown environment variable!")
    else
      select case (ctype)
      case (COLLCT_INTEGER)
        call output_line (trim(sys_siL(&
          collct_getvalue_int (rcmdStatus%rcollection, svalname) ,10)))
          
      case (COLLCT_REAL)
        call output_line (trim(sys_sdEL(&
          collct_getvalue_real (rcmdStatus%rcollection, svalname) ,10)))
          
      case (COLLCT_STRING)
        call collct_getvalue_string (rcmdStatus%rcollection, svalname, svalue)
        call output_line (svalue)
          
      case default
        call output_line ("No detailed information available.")
      end select
      
    end if
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_info (rcmdStatus,Sargs)
  
  !<description>
    ! Command: INFO.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: svalname
    character(len=SYS_STRLEN) :: svalue
    integer :: ctype
    logical :: bexists
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the next token(s)
    call cmdprs_getparam (Sargs,2,svalname,.true.,.false.)
    
    ! That's an awful task. Figure out the type at first.
    ! Then get the value and/or print information.
    ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists)
    
    if (.not. bexists) then
      call output_line ("Unknown environment variable!")
    else
      call output_line ("Variable: "//svalname)
      select case (ctype)
      case (COLLCT_INTEGER)
        call output_line ("Type    : Integer")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call output_line (trim(sys_siL(&
          collct_getvalue_int (rcmdStatus%rcollection, svalname) ,10)))
          
      case (COLLCT_REAL)
        call output_line ("Type    : Double precision")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call output_line (trim(sys_sdEL(&
          collct_getvalue_real (rcmdStatus%rcollection, svalname) ,10)))
          
      case (COLLCT_STRING)
        call output_line ("Type    : string")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call collct_getvalue_string (rcmdStatus%rcollection, svalname, svalue)
        call output_line (svalue)
          
      case (COLLCT_BOUNDARY)
        call output_line ("Type    : Boundary object (2D)")

      case (COLLCT_TRIA)
        call output_line ("Type    : Triangulation object")
        call output_lbrk ()
        p_rtriangulation => collct_getvalue_tria (rcmdStatus%rcollection, svalname)
        call tria_infoStatistics (p_rtriangulation,.true.)

      case (COLLCT_MSHHIERARCHY)
        call output_line ("Type    : Mesh hierarchy")
        call output_lbrk ()
        p_rmeshhierarchy => collct_getvalue_mshh (rcmdStatus%rcollection, svalname)
        call mshh_printHierStatistics (p_rmeshhierarchy)

      case (COLLCT_FESPACE)
        call output_line ("Type    : FE space")
        call output_lbrk ()
        p_rfeSpace => collct_getvalue_fesp (rcmdStatus%rcollection, svalname)
        call fesph_infoStatistics (p_rfeSpace,.true.)
        call spdiscr_infoBlockDiscr (p_rfeSpace%p_rdiscretisation)

      case (COLLCT_FEHIERARCHY)
        call output_line ("Type    : Mesh hierarchy")
        call output_lbrk ()
        p_rfeHierarchy => collct_getvalue_feh (rcmdStatus%rcollection, svalname)
        call fesph_printHierStatistics (p_rfeHierarchy)

      case (COLLCT_MLPRJHIERARCHY)
        call output_line ("Type    : Multilevel projection hierarchy")
        call output_lbrk ()

      case (COLLCT_BLKVECTOR)
        call output_line ("Type    : Block vector")
        call output_lbrk ()
        p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svalname)
        call lsysbl_infoVector (p_rvectorBlock)

      case default
        call output_line ("Type    : "//sys_siL(ctype,10))
      end select
      
    end if
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_read2dprm (rcmdStatus,Sargs)
  
  !<description>
    ! Command: READ2DPRM.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sfilename
    character(len=COLLCT_MLNAME) :: sname
    type(t_boundary), pointer :: p_rboundary
  
    if (size(Sargs) .lt. 3) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the identifier and file name 
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    call cmdprs_getparam (Sargs,3,sfilename,.false.,.true.)

    ! Open the file and read.
    inquire(file=trim(sfilename), exist=bexists)
    if (.not. bexists) then
      call output_line ("File not found!")
    else
      call output_line ("Reading file: "//trim(sfilename))
    
      ! Read
      allocate (p_rboundary)
      call boundary_read_prm(p_rboundary, sfilename)
      
      ! Remove old value from collection if present
      call do_destroy(rcmdStatus,sname,.false.)
      
      ! Add to the collection
      call collct_setvalue_bdry (rcmdStatus%rcollection, sname, p_rboundary, .true.)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_read2dtri (rcmdStatus,Sargs)
  
  !<description>
    ! Command: READ2DTRI.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sfilename, stoken
    character(len=COLLCT_MLNAME) :: sname
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
  
    if (size(Sargs) .lt. 3) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the identifier and file name 
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    call cmdprs_getparam (Sargs,3,sfilename,.false.,.true.)

    ! Open the file and read.
    inquire(file=trim(sfilename), exist=bexists)
    if (.not. bexists) then
      call output_line ("File not found!")
    else
      nullify(p_rboundary)
    
      ! Is there a boundary object attached?
      if (size(Sargs) .ge. 5) then
        call cmdprs_getparam (Sargs,4,stoken,.true.,.false.)
        if (stoken .eq. "--BOUNDARY") then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,5,stoken,.false.,.true.)
          
          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        end if
      end if
    
      call output_line ("Reading file: "//trim(sfilename))
    
      ! Read
      allocate (p_rtriangulation)
      
      if (associated(p_rboundary)) then
        call tria_readTriFile2D(p_rtriangulation, sfilename, p_rboundary)
        
        ! Create a standard mesh.
        call tria_initStandardMeshFromRaw (p_rtriangulation, p_rboundary)
      else
        call output_line("Warning. No boundary specified!")
      
        call tria_readTriFile2D(p_rtriangulation, sfilename)
        
        ! Create a standard mesh.
        call tria_initStandardMeshFromRaw (p_rtriangulation)
      end if
      
      ! Remove old value from collection if present
      call do_destroy(rcmdStatus,sname,.false.)
      
      ! Add to the collection
      call collct_setvalue_tria (rcmdStatus%rcollection, sname, p_rtriangulation, .true.)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_meshrefine (rcmdStatus,Sargs)
  
  !<description>
    ! Command: MESHREFINE.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: stoken
    integer :: cmethod, icount
    integer :: i
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the identifier and file name 
    call cmdprs_getparam (Sargs,2,stoken,.true.,.false.)
    p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,bexists=bexists)
    if (.not. bexists) then
      call output_line ("Unknown variable!")
      return
    end if
    
    ! Read the other parameters
    cmethod = 0 ! 2-level ordering
    icount = 1
    nullify(p_rboundary)
  
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--BOUNDARY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--TIMES") then
        if (i .le. size(Sargs)) then
          ! Get the method type
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          ! Number of refinements
          read (stoken,*) icount
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--METHOD") then
        if (i .le. size(Sargs)) then
          ! Get the method type
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          if (stoken .eq. "--2LEVELORDERED") then
            ! 2-level ordering
            cmethod = 0
          else
            call output_line ("Warning: Unknown refinement method! Using default.")
            cmethod = 0
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      end if
    end do
  
    call output_line ("Refining mesh... [",bnolinebreak=.true.)
  
    ! Refine the mesh
    select case (cmethod)
    case (0)
      ! 2-level ordered
      
      do i = 1,icount
        call output_line (" "//trim(sys_siL(i,10)),bnolinebreak=.true.)

        ! Boundary present?
        if (associated(p_rboundary)) then
          call tria_quickRefine2LevelOrdering(1,p_rtriangulation)
        else
          call tria_quickRefine2LevelOrdering(1,p_rtriangulation,p_rboundary)
        end if
      end do

      ! Create a standard mesh.
      if (associated(p_rboundary)) then
        call tria_initStandardMeshFromRaw (p_rtriangulation)
      else
        call tria_initStandardMeshFromRaw (p_rtriangulation,p_rboundary)
      end if
      
      call output_line ("]")
    end select
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_meshhierarchy (rcmdStatus,Sargs)
  
  !<description>
    ! Command: MESHHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken
    integer :: cmethod, ilevels
    integer :: i
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_meshHierarchy), pointer :: p_rmeshHierarchy
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    cmethod = 0 ! 2-level ordering
    ilevels = 1
    nullify(p_rboundary)
    nullify(p_rtriangulation)
  
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--BOUNDARY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--MESH") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Mesh object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--LEVELS") then
        if (i .le. size(Sargs)) then
          ! Get the method type
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          ! Number of refinements
          read (stoken,*) ilevels
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--METHOD") then
        if (i .le. size(Sargs)) then
          ! Get the method type
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          if (stoken .eq. "2LEVELORDERED") then
            ! 2-level ordering
            cmethod = 0
          else
            call output_line ("Warning: Unknown refinement method! Using default.")
            cmethod = 0
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      end if
    end do
    
    if (.not. associated (p_rtriangulation)) then
      call output_line ("Invalid triangulation!")
      return
    end if
  
    ! Create the hierarchy.
    allocate(p_rmeshHierarchy)
    
    ! Create the hierarchy
    select case (cmethod)
    case (0)
      ! 2-level ordered
  
      ! Boundary present?
      if (associated(p_rboundary)) then
        call output_line ("Creating mesh hierarchy... [1",bnolinebreak=.true.)
        call mshh_initHierarchy (p_rmeshHierarchy,p_rtriangulation,0,ilevels,&
            rboundary=p_rboundary)
            
        call mshh_refineHierarchy2lv (p_rmeshHierarchy,ilevels,&
            rboundary=p_rboundary,bprint=.true.)
      else
        call output_line ("Warning: No boundary present!")
        
        call output_line ("Creating mesh hierarchy... [1",bnolinebreak=.true.)
        call mshh_initHierarchy (p_rmeshHierarchy,p_rtriangulation,0,ilevels)

        call mshh_refineHierarchy2lv (p_rmeshHierarchy,ilevels,bprint=.true.)
      end if

      call output_line ("]")
    end select
    
    ! Remove old value from collection if present
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_mshh (rcmdStatus%rcollection, sname, p_rmeshHierarchy, .true.)
    
  end subroutine

  ! ***************************************************************************

    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    
    !<description>
      ! Callback routine to set up a block discretisation.
    !</description>
    
    !<input>
      ! Triangulation structure
      type(t_triangulation), intent(in) :: rtriangulation
  
      ! Definition of the boundary
      type(t_boundary), intent(in), optional :: rboundary
    
      ! Collection structure with information about the discretisation
      type(t_collection), intent(inout), optional :: rcollection
    !</input>
  
    !<output>
      ! Block discretisation structure
      type(t_blockDiscretisation), intent(out) :: rdiscr
    !</output>

      ! local variables
      integer :: i,nspaces
      integer, dimension(:), allocatable :: p_IelementIdsTri
      integer, dimension(:), allocatable :: p_IelementIdsQuad
      integer, dimension(:), allocatable :: p_IcubIdsTri
      integer, dimension(:), allocatable :: p_IcubIdsQuad

      nspaces = rcollection%IquickAccess(1)
      allocate(p_IelementIdsTri(nspaces))
      allocate(p_IelementIdsQuad(nspaces))
      allocate(p_IcubIdsTri(nspaces))
      allocate(p_IcubIdsQuad(nspaces))
      
      call collct_getvalue_intarr (rcollection, "IelementIdsTri", p_IelementIdsTri)
      call collct_getvalue_intarr (rcollection, "IelementIdsQuad", p_IelementIdsQuad)
      call collct_getvalue_intarr (rcollection, "IcubIdsTri", p_IcubIdsTri)
      call collct_getvalue_intarr (rcollection, "IcubIdsQuad", p_IcubIdsQuad)

      select case (rtriangulation%ndim)
      case (NDIM2D)
        ! Create a block discretisation of the specific size.
        call spdiscr_initBlockDiscr (rdiscr,nspaces,rtriangulation,rboundary)
        
        ! Create the sub-discretisation structures.
        do i=1,nspaces
          if (p_IelementIdsQuad(i) .eq. 0) then
            ! Pure tri space
            call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(i), &
                INT(p_IelementIdsTri(i),I32),INT(p_IcubIdsTri(i),I32),rtriangulation, &
                rboundary)
                
          else if (p_IelementIdsTri(i) .eq. 0) then
            ! Pure quad space
            call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(i), &
                INT(p_IelementIdsQuad(i),I32),INT(p_IcubIdsQuad(i),I32),rtriangulation, &
                rboundary)
          else
            ! Combined space
            call spdiscr_initDiscr_triquad (rdiscr%RspatialDiscr(i), &
                int(p_IelementIdsTri(i),I32), INT(p_IelementIdsQuad(i),I32),&
                int(p_IcubIdsTri(i),I32), INT(p_IcubIdsQuad(i),I32),rtriangulation, &
                rboundary)
          end if
        end do
        
      case default
        call output_line("Error. Dimension not supported!")
        call sys_halt()
      end select

    end subroutine

  ! ***************************************************************************

    subroutine prepare_fgetDiscr(rcollection,rcmdStatus,nspaces,&
        IelementIdsTri,IelementIdsQuad,IcubIdsTri,IcubIdsQuad)
    
    !<description>
      ! Prepares setting up a discretisation.
      ! Must be called in advance to fgetDiscr.
    !</description>
    
    !<input>
      ! Number of components
      integer, intent(in) :: nspaces
    
      ! List of element id's for tri/tetra elements
      integer, dimension(:), intent(in), target :: IelementIdsTri

      ! List of element id's for quad/hexa elements
      integer, dimension(:), intent(in), target :: IelementIdsQuad

      ! List of cubature id's for tri/tetra elements
      integer, dimension(:), intent(in), target :: IcubIdsTri

      ! List of cubature id's for quad/hexa elements
      integer, dimension(:), intent(in), target :: IcubIdsQuad
    !</input>
    
    !<inputoutput>
      ! Current status block.
      type(t_commandstatus), intent(inout), target :: rcmdStatus

      ! Collection structure with information about the discretisation
      type(t_collection), intent(inout), optional :: rcollection
    !</inputoutput>

      rcollection%p_rnextCollection => rcmdStatus%rcollection
      
      rcollection%IquickAccess(1) = nspaces
      call collct_setvalue_intarr (rcollection, "IelementIdsTri", IelementIdsTri, .true.)
      call collct_setvalue_intarr (rcollection, "IelementIdsQuad", IelementIdsQuad, .true.)
      call collct_setvalue_intarr (rcollection, "IcubIdsTri", IcubIdsTri, .true.)
      call collct_setvalue_intarr (rcollection, "IcubIdsQuad", IcubIdsQuad, .true.)
      
    end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_fespace (rcmdStatus,Sargs)
  
  !<description>
    ! Command: FESPACE.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken
    integer :: cmethod, ilevels, ncomponents
    integer :: i,j
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_feSpaceLevel), pointer :: p_rfeSpace1,p_rfeSpace2
    type (t_collection) :: rcollection
    integer, dimension(:), allocatable :: p_IelementIdsTri,p_IelementIdsQuad
    integer, dimension(:), allocatable :: p_IcubIdsTri,p_IcubIdsQuad
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    cmethod = 0 ! 2-level ordering
    ncomponents = 0
    ilevels = 1
    nullify(p_rboundary)
    nullify(p_rtriangulation)
    
    ! For concatenation
    nullify(p_rfeSpace1)
    nullify(p_rfeSpace2)
  
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--BOUNDARY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--MESH") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Mesh object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--COMPONENTS") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          read (stoken,*) ncomponents
          
          if (ncomponents .le. 0) then
            call output_line("Error. Number of components invalid!")
            return
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--TRIELEMENTS") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IelementIdsTri(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IelementIdsTri(j) = elem_igetID(stoken,.true.)
            if (p_IelementIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              return
            end if

            if (elem_igetDimension(p_IelementIdsTri(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              return
            end if

            if (elem_igetNVE(p_IelementIdsTri(j)) .ne. 3) then
              call output_line("Error. Not a tri element: "//trim(stoken))
              return
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--QUADELEMENTS") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IelementIdsQuad(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IelementIdsQuad(j) = elem_igetID(stoken,.true.)
            if (p_IelementIdsQuad(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              return
            end if

            if (elem_igetDimension(p_IelementIdsQuad(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              return
            end if

            if (elem_igetNVE(p_IelementIdsQuad(j)) .ne. 4) then
              call output_line("Error. Not a quad element: "//trim(stoken))
              return
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--TRICUB") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IcubIdsTri(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IcubIdsTri(j) = cub_igetID(stoken,.true.)
            if (p_IelementIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              return
            end if

          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--QUADCUB") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IcubIdsQuad(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IcubIdsQuad(j) = cub_igetID(stoken,.true.)
            if (p_IelementIdsQuad(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              return
            end if

          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--CONCAT") then
        if (i .lt. size(Sargs)-1) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          if (stoken .eq. sname) then
            call output_line("Error. Destination variable must be different from source!")
            return
          end if

          p_rfeSpace1 => collct_getvalue_fesp(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            call output_line("Error. 1st source FE space does not exist.")
            return
          end if

          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          if (stoken .eq. sname) then
            call output_line("Error. Destination variable must be different from source!")
            return
          end if

          p_rfeSpace2 => collct_getvalue_fesp(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            call output_line("Error. 2nd source FE space does not exist.")
            return
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rtriangulation)) then
      call output_line ("Invalid triangulation!")
      return
    end if
  
    ! Create the fe space.
    allocate(p_rfeSpace)
    
    if (associated (p_rfeSpace1)) then
      ! Concat two predefined FE spaces.
      call fesph_concatFeSpaces (p_rfeSpace1,p_rfeSpace2,p_rfeSpace,p_rtriangulation)
    else
      if ((.not. allocated(p_IelementIdsTri)) .and. &
          (.not. allocated(p_IelementIdsQuad))) then
        call output_line("Error. No element ID's specified!")
        return
      end if
    
      ! Create missing id arrays.
      if (.not. allocated(p_IelementIdsTri)) then
        allocate(p_IelementIdsTri(ncomponents))
        p_IelementIdsTri(:) = 0
      end if
      
      if (.not. allocated(p_IelementIdsQuad)) then
        allocate(p_IelementIdsQuad(ncomponents))
        p_IelementIdsQuad(:) = 0
      end if
    
      if (.not. allocated(p_IcubIdsTri)) then
        allocate(p_IcubIdsTri(ncomponents))
        p_IcubIdsTri(:) = 0
      end if

      ! Probably take the standard cubature formula.      
      where ((p_IcubIdsTri .eq. 0) .and. (p_IelementIdsTri .ne. 0)) 
        p_IcubIdsTri = spdiscr_getStdCubature(int(p_IelementIdsTri,I32))
      end where

      ! Probably take the standard cubature formula.      
      where ((p_IcubIdsQuad .eq. 0) .and. (p_IelementIdsQuad .ne. 0)) 
        p_IcubIdsQuad = spdiscr_getStdCubature(int(p_IelementIdsQuad,I32))
      end where
    
      if (.not. allocated(p_IcubIdsQuad)) then
        allocate(p_IcubIdsQuad(ncomponents))
        p_IcubIdsQuad(:) = 0
      end if
    
      ! Create the FE space usinf fgetDiscr.
      call collct_init (rcollection)
      call prepare_fgetDiscr(rcollection,rcmdStatus,ncomponents,&
          p_IelementIdsTri,p_IelementIdsQuad,p_IcubIdsTri,p_IcubIdsQuad)
      if (associated(p_rboundary)) then
        call fesph_createFEspace (p_rfeSpace,1,&
            p_rtriangulation,1,fgetDiscr,rcollection,rboundary=p_rboundary)
      else
        call output_line ("Warning: No boundary present!")
        
        call fesph_createFEspace (p_rfeSpace,1,&
            p_rtriangulation,1,fgetDiscr,rcollection)
      end if
      call collct_done (rcollection)
    end if
    
    ! Remove old value from collection if present
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_fesp (rcmdStatus%rcollection, sname, p_rfeSpace, .true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_fehierarchy (rcmdStatus,Sargs)
  
  !<description>
    ! Command: FEHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken
    integer :: cmethod, ilevels, ncomponents
    integer :: i,j
    type(t_boundary), pointer :: p_rboundary
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type (t_collection) :: rcollection
    integer, dimension(:), allocatable :: p_IelementIdsTri,p_IelementIdsQuad
    integer, dimension(:), allocatable :: p_IcubIdsTri,p_IcubIdsQuad
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    cmethod = 0 ! 2-level ordering
    ncomponents = 0
    ilevels = 1
    nullify(p_rboundary)
    nullify(p_rmeshhierarchy)
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--BOUNDARY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--MESHHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rmeshhierarchy => collct_getvalue_mshh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Mesh hierarchy object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--COMPONENTS") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          read (stoken,*) ncomponents
          
          if (ncomponents .le. 0) then
            call output_line("Error. Number of components invalid!")
            return
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--TRIELEMENTS") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IelementIdsTri(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IelementIdsTri(j) = elem_igetID(stoken,.true.)
            if (p_IelementIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              return
            end if

            if (elem_igetDimension(p_IelementIdsTri(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              return
            end if

            if (elem_igetNVE(p_IelementIdsTri(j)) .ne. 3) then
              call output_line("Error. Not a tri element: "//trim(stoken))
              return
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--QUADELEMENTS") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IelementIdsQuad(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IelementIdsQuad(j) = elem_igetID(stoken,.true.)
            if (p_IelementIdsQuad(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              return
            end if

            if (elem_igetDimension(p_IelementIdsQuad(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              return
            end if

            if (elem_igetNVE(p_IelementIdsQuad(j)) .ne. 4) then
              call output_line("Error. Not a quad element: "//trim(stoken))
              return
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--TRICUB") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IcubIdsTri(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IcubIdsTri(j) = cub_igetID(stoken,.true.)
            if (p_IelementIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              return
            end if

          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--QUADCUB") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IcubIdsQuad(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IcubIdsQuad(j) = cub_igetID(stoken,.true.)
            if (p_IelementIdsQuad(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              return
            end if

          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rmeshhierarchy)) then
      call output_line ("Invalid mesh hierarchy!")
      return
    end if
  
    ! Create the fe space.
    allocate(p_rfeHierarchy)
    
    if ((.not. allocated(p_IelementIdsTri)) .and. &
        (.not. allocated(p_IelementIdsQuad))) then
      call output_line("Error. No element ID's specified!")
      return
    end if
  
    ! Create missing id arrays.
    if (.not. allocated(p_IelementIdsTri)) then
      allocate(p_IelementIdsTri(ncomponents))
      p_IelementIdsTri(:) = 0
    end if
    
    if (.not. allocated(p_IelementIdsQuad)) then
      allocate(p_IelementIdsQuad(ncomponents))
      p_IelementIdsQuad(:) = 0
    end if
  
    if (.not. allocated(p_IcubIdsTri)) then
      allocate(p_IcubIdsTri(ncomponents))
      p_IcubIdsTri(:) = 0
    end if

    ! Probably take the standard cubature formula.      
    where ((p_IcubIdsTri .eq. 0) .and. (p_IelementIdsTri .ne. 0)) 
      p_IcubIdsTri = spdiscr_getStdCubature(int(p_IelementIdsTri,I32))
    end where

    ! Probably take the standard cubature formula.      
    where ((p_IcubIdsQuad .eq. 0) .and. (p_IelementIdsQuad .ne. 0)) 
      p_IcubIdsQuad = spdiscr_getStdCubature(int(p_IelementIdsQuad,I32))
    end where
  
    if (.not. allocated(p_IcubIdsQuad)) then
      allocate(p_IcubIdsQuad(ncomponents))
      p_IcubIdsQuad(:) = 0
    end if
  
    call output_line ("Creating FE hierarchy.")
  
    ! Create the FE space usinf fgetDiscr.
    call collct_init (rcollection)
    call prepare_fgetDiscr(rcollection,rcmdStatus,ncomponents,&
        p_IelementIdsTri,p_IelementIdsQuad,p_IcubIdsTri,p_IcubIdsQuad)
    if (associated(p_rboundary)) then
      call fesph_createHierarchy (p_rfeHierarchy,p_rmeshHierarchy%nlevels,&
          p_rmeshHierarchy,fgetDiscr,rcollection,rboundary=p_rboundary)
    else
      call output_line ("Warning: No boundary present!")
      
      call fesph_createHierarchy (p_rfeHierarchy,p_rmeshHierarchy%nlevels,&
          p_rmeshHierarchy,fgetDiscr,rcollection)
    end if
    call collct_done (rcollection)
    
    ! Remove old value from collection if present
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_feh (rcmdStatus%rcollection, sname, p_rfeHierarchy, .true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_mlevelprjhierarchy (rcmdStatus,Sargs)
  
  !<description>
    ! Command: MLEVELPRJHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken
    integer :: i
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_interlevelProjectionHier), pointer :: p_rprjHier
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    nullify(p_rfeHierarchy)
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--FEHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rfeHierarchy => collct_getvalue_feh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rfeHierarchy)
            call output_line("Warning. FE hierarchy object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rfeHierarchy)) then
      call output_line ("Invalid mesh hierarchy!")
      return
    end if
  
    ! Create the fe space.
    allocate(p_rprjHier)
    
    call output_line ("Creating multilevel projection hierarchy.")
  
    ! Create the FE space usinf fgetDiscr.
    call mlprj_initPrjHierarchy(p_rprjHier,1,p_rfeHierarchy%nlevels)
    do i = 1,p_rfeHierarchy%nlevels
      call mlprj_initPrjHierarchyLevel(p_rprjHier,i,&
          p_rfeHierarchy%p_rfeSpaces(i)%p_rdiscretisation)
    end do
    call mlprj_commitPrjHierarchy (p_rprjHier)

    ! Remove old value from collection if present
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_mlprjh (rcmdStatus%rcollection, sname, p_rprjHier, .true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_createblockvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: READBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken, sfilename
    integer :: i,ilevel
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
    logical :: bformatted
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    ilevel = 0
    nullify(p_rfeHierarchy)
    bformatted = .true.
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--FEHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          p_rfeHierarchy => collct_getvalue_feh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rfeHierarchy)
            call output_line("Warning. FE hierarchy object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--LEVEL") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          read (stoken,*) ilevel
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rfeHierarchy)) then
      call output_line ("Invalid mesh hierarchy!")
      return
    end if
    
    if (ilevel .eq. 0) ilevel = p_rfeHierarchy%nlevels
  
    inquire(file=trim(sfilename), exist=bexists)
    
    call output_line ("Creating block vector.")

    ! Create the vector.
    allocate(p_rvectorBlock)
    call lsysbl_createVectorBlock (p_rfeHierarchy%p_rfeSpaces(ilevel)%p_rdiscretisation,&
        p_rvectorBlock,.true.)

    ! Remove old value from collection if present
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_vec (rcmdStatus%rcollection, sname, p_rvectorBlock, .true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_readblockvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: READBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken, sfilename
    integer :: i,ilevel
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
    logical :: bformatted
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    ilevel = 1
    nullify(p_rfeHierarchy)
    bformatted = .true.
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--FILENAME") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          sfilename = trim(stoken)
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--UNFORMATTED") then

        bformatted = .false.

      end if      

    end do
    
    if (sfilename .eq. "") then
      call output_line ("Invalid filename!")
      return
    end if

    inquire(file=trim(sfilename), exist=bexists)
    
    if (.not. bexists) then
      call output_line ("File not found!")
    end if
    
    p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, sname)
    if (.not. associated(p_rvectorBlock)) then
      call output_line ("Invalid vector!")
      return
    end if
    
    call output_line ("Reading vector: "//trim(sfilename))

    ! Create the vector.
    call vecio_readBlockVectorHR (p_rvectorBlock, stoken, .true.,&
        0, sfilename, bformatted)

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_writeblockvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: WRITEBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=SYS_STRLEN) :: sname, stoken, sfilename, sformat
    integer :: i,ilevel
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
    logical :: bformatted
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    ilevel = 1
    nullify(p_rfeHierarchy)
    bformatted = .true.
    sformat = "(E20.10)"
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--FILENAME") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          sfilename = trim(stoken)
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--UNFORMATTED") then

        bformatted = .false.

      else if (stoken .eq. "--FORMAT") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          sformat = "("//trim(adjustl(stoken))//")"
          
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (sfilename .eq. "") then
      call output_line ("Invalid filename!")
      return
    end if
    
    p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, sname)
    if (.not. associated(p_rvectorBlock)) then
      call output_line ("Invalid vector!")
      return
    end if
    
    call output_line ("Writing vector: "//trim(sfilename))
    
    if (bformatted) then
      call vecio_writeBlockVectorHR (p_rvectorBlock, "vector", .true.,&
          0, sfilename, sformat)
    else
      call vecio_writeBlockVectorHR (p_rvectorBlock, "vector", .true.,&
          0, sfilename)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_copyvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: WRITEBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    character(len=SYS_STRLEN) :: ssourceall, sdestall, ssource, sdest
    integer, dimension(:), pointer :: p_IdxSource,p_IdxDest
    logical :: bsourcescalar, bdestscalar
    integer :: icomp1, icomp2
    
    if (size(Sargs) .lt. 3) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Source and destination vector
    call cmdprs_getparam (Sargs,2,ssourceall,.true.,.false.)
    call cmdprs_getparam (Sargs,3,sdestall,.true.,.false.)
    
    ! Vector name and components
    call cmdprs_parseIndexString1D (ssourceall,ssource,p_IdxSource)
    call cmdprs_parseIndexString1D (sdestall,sdest,p_IdxDest)
    
    bsourcescalar = associated(p_IdxSource)
    bdestscalar = associated(p_IdxDest)
    if (bsourcescalar .neqv. bdestscalar) then
      if (bsourcescalar) deallocate(p_IdxSource)
      if (bdestscalar) deallocate(p_IdxDest)
      call output_line ("Both vectors must be block or scalar")
      return
    end if

    p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, ssource)
    if (.not. associated(p_rvectorBlock1)) then
      if (bsourcescalar) deallocate(p_IdxSource)
      if (bdestscalar) deallocate(p_IdxDest)
      call output_line ("Invalid source vector!")
      return
    end if
    
    p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, sdest)
    if (.not. associated(p_rvectorBlock2)) then
      if (bsourcescalar) deallocate(p_IdxSource)
      if (bdestscalar) deallocate(p_IdxDest)
      call output_line ("Invalid destination vector!")
      return
    end if

    if (.not. bsourcescalar) then
      ! Block copy
      call lsysbl_copyVector (p_rvectorBlock1,p_rvectorBlock2)
    else
      if ((ubound(p_IdxSource,1) .ne. 1) .or. (ubound(p_IdxDest,1) .ne. 1)) then
        deallocate(p_IdxSource)
        deallocate(p_IdxDest)
        call output_line ("Only one component allowed!")
        return
      end if
      
      icomp1 = p_IdxSource(1)
      icomp2 = p_IdxDest(1)
    
      if ((icomp1 .lt. 1) .or. (icomp1 .gt. p_rvectorBlock1%nblocks)) then
        deallocate(p_IdxSource)
        deallocate(p_IdxDest)
        call output_line ("Invalid source component!")
        return
      end if

      if ((icomp2 .lt. 1) .or. (icomp2 .gt. p_rvectorBlock1%nblocks)) then
        deallocate(p_IdxSource)
        deallocate(p_IdxDest)
        call output_line ("Invalid destination component!")
        return
      end if
      
      call lsyssc_copyVector (p_rvectorBlock1%RvectorBlock(icomp1),&
          p_rvectorBlock2%RvectorBlock(icomp2))

      deallocate(p_IdxSource)
      deallocate(p_IdxDest)
    end if
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_interpolatevector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: INTERPOLATEVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    character(len=COLLCT_MLNAME) :: ssource, sdest, stoken
    type(t_interlevelProjectionHier), pointer :: p_rprjHier
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    integer :: i,isource,idest
    type(t_vectorBlock), pointer :: p_rtemp1,p_rtemp2
    
    if (size(Sargs) .lt. 6) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Source vector and component
    call cmdprs_getparam (Sargs,2,ssource,.true.,.false.)
    call cmdprs_getparam (Sargs,3,stoken,.true.,.false.)
    read (stoken,*) isource
    
    ! Destination
    call cmdprs_getparam (Sargs,4,sdest,.true.,.false.)
    call cmdprs_getparam (Sargs,5,stoken,.true.,.false.)
    read (stoken,*) idest

    p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, ssource)
    if (.not. associated(p_rvectorBlock1)) then
      call output_line ("Invalid source vector!")
      return
    end if

    p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, sdest)
    if (.not. associated(p_rvectorBlock2)) then
      call output_line ("Invalid destination vector!")
      return
    end if

    nullify(p_rfeHierarchy)
    nullify(p_rprjHier)

    i = 5
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--MLPRJHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rprjHier => collct_getvalue_mlprjh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rprjHier)
            call output_line("Warning. Multilevel projection hierarchy does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--FEHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          p_rfeHierarchy => collct_getvalue_feh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rfeHierarchy)
            call output_line("Warning. FE hierarchy object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rprjHier)) then
      call output_line ("Invalid multilevel projection hierarchy!")
      return
    end if
    
    if ((isource .lt. p_rprjHier%nlmin) .or. (isource .gt. p_rprjHier%nlmax)) then
      call output_line ("Invalid source level!")
      return
    end if
    
    if ((idest .lt. p_rprjHier%nlmin) .or. (idest .gt. p_rprjHier%nlmax)) then
      call output_line ("Invalid destination level!")
      return
    end if
    
    ! Go either up or down - or copy.
    if (isource .eq. idest) then
    
      call lsysbl_copyVector (p_rvectorBlock1,p_rvectorBlock2)
      
    else if (isource .lt. idest) then
    
      if (isource .eq. idest-1) then
      
        ! One prolongation
        p_rtemp1 => p_rvectorBlock1
      
      else
      
        if (.not. associated (p_rfeHierarchy)) then
          call output_line ("Invalid FE hierarchy!")
          return
        end if
      
        ! Multiple prolongations
        allocate (p_rtemp1)
        call lsysbl_createVectorBlock (p_rvectorBlock1,p_rtemp1,.false.)
        call lsysbl_copyVector (p_rvectorBlock1,p_rtemp1)
        
        ! Create the temp vectors using the FE hierarchy.
        do i=isource+1,idest-1
          allocate (p_rtemp2)
          call lsysbl_createVectorBlock (p_rfeHierarchy%p_rfeSpaces(i)%p_rdiscretisation,&
              p_rtemp2,.false.)
          call mlprj_performProlongationHier (p_rprjHier,&
              i,p_rtemp1,p_rtemp2)
          call lsysbl_releaseVector (p_rtemp1)
          deallocate (p_rtemp1)
          p_rtemp1 => p_rtemp2
        end do
      end if

      ! Final prolongation
      call mlprj_performProlongationHier (p_rprjHier,&
          idest,p_rtemp1,p_rvectorBlock2)
          
      if (isource .eq. idest-1) then
        ! Cleanup
        call lsysbl_releaseVector (p_rtemp1)
        deallocate(p_rtemp1)
      end if
    
    else
      ! Interpolation. NOT RESTRICTION!!!
    
      if (isource-1 .eq. idest) then
      
        ! One prolongation
        p_rtemp1 => p_rvectorBlock1
      
      else
      
        if (.not. associated (p_rfeHierarchy)) then
          call output_line ("Invalid FE hierarchy!")
          return
        end if
        
        ! Multiple interpolations
        allocate (p_rtemp1)
        call lsysbl_createVectorBlock (p_rvectorBlock1,p_rtemp1,.false.)
        call lsysbl_copyVector (p_rvectorBlock1,p_rtemp1)
        
        ! Create the temp vectors using the FE hierarchy.
        do i=idest-1,isource+1,-1
          allocate (p_rtemp2)
          call lsysbl_createVectorBlock (p_rfeHierarchy%p_rfeSpaces(i)%p_rdiscretisation,&
              p_rtemp2,.false.)
          call mlprj_performInterpolationHier (p_rprjHier,&
              i+1,p_rtemp2,p_rtemp1)
          call lsysbl_releaseVector (p_rtemp1)
          deallocate (p_rtemp1)
          p_rtemp1 => p_rtemp2
        end do
      end if

      ! Final interpolation
      call mlprj_performInterpolationHier (p_rprjHier,&
          idest+1,p_rvectorBlock2,p_rtemp1)
          
      if (isource-1 .eq. idest) then
        ! Cleanup
        call lsysbl_releaseVector (p_rtemp1)
        deallocate(p_rtemp1)
      end if
      
    end if

  end subroutine

  ! ***************************************************************************

    subroutine fcoeff_analytPrj (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, &
                  IdofsTest, rdomainIntSubset, &
                  Dcoefficients, rcollection)
    
    ! Returns values of an FE function in cubature points.
    
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    type(t_linearForm), intent(in) :: rform
    integer, intent(in) :: nelements
    integer, intent(in) :: npointsPerElement
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    integer, dimension(:,:), intent(in) :: IdofsTest
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset    
    type(t_collection), intent(inout), optional :: rcollection
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

      ! local variables
      integer :: icomponent
      type(t_vectorBlock), pointer :: p_rvectorBlock
      
      ! Get the component and the FE function
      p_rvectorBlock => rcollection%p_rvectorQuickAccess1
      icomponent = rcollection%IquickAccess(1)

      ! Evaluate the FE function      
      call fevl_evaluate_sim (p_rvectorBlock%RvectorBlock(icomponent), &
          rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
  
    end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_l2projection (rcmdStatus,Sargs)
  
  !<description>
    ! Command: INTERPOLATEVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    character(len=COLLCT_MLNAME) :: ssource, sdest, stoken
    type(t_matrixScalar) :: rmatrixMass
    integer :: i
    logical :: bverbose
    type(t_collection) :: rcollection
    type(t_configL2ProjectionByMass) :: rL2ProjectionConfig
    
    if (size(Sargs) .lt. 3) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Source vector and destination
    call cmdprs_getparam (Sargs,2,ssource,.true.,.false.)
    call cmdprs_getparam (Sargs,3,sdest,.true.,.false.)

    p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, ssource)
    if (.not. associated(p_rvectorBlock1)) then
      call output_line ("Invalid source vector!")
      return
    end if

    p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, sdest)
    if (.not. associated(p_rvectorBlock2)) then
      call output_line ("Invalid destination vector!")
      return
    end if

    i = 4
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--VERBOSE") then
        bverbose = .true.
      end if      

      if (stoken .eq. "--RELERROR") then
        if (i .le. size(Sargs)) then
          call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
          i = i+1
          read (stoken,*) rL2ProjectionConfig%depsrel
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      end if      

    end do

    ! Clear the destination
    call lsysbl_clearVector (p_rvectorBlock2)

    ! Loop through the components
    do i=1,min(p_rvectorBlock1%nblocks,p_rvectorBlock2%nblocks)
      if (bverbose) then
        call output_lbrk ()
        call output_line ("Component : "//trim(sys_siL(i,10)))
        call output_line ("Creating mass matrix - structure...")
      end if

      ! Create a mass matrix in that space
      call bilf_createMatrixStructure (p_rvectorBlock2%p_rblockDiscr%RspatialDiscr(i),&
          LSYSSC_MATRIX9,rmatrixMass)

      if (bverbose) then
        call output_line ("Creating mass matrix - content...")
      end if

      call stdop_assembleSimpleMatrix (rmatrixMass,DER_FUNC,DER_FUNC,1.0_DP,.true.)
      
      if (bverbose) then
        call output_line ("Projecting...")
      end if

      ! Do the L2 projection
      rcollection%p_rvectorQuickAccess1 => p_rvectorBlock1
      rcollection%IquickAccess(1) = i
      call anprj_analytL2projectionByMass (p_rvectorBlock2%RvectorBlock(i), rmatrixMass,&
          fcoeff_analytPrj, rcollection, rL2ProjectionConfig)

      if (bverbose) then
        call output_line ("Rel. error: "//trim(sys_sdEL(rL2ProjectionConfig%drelError,10)))
        call output_line ("Abs. error: "//trim(sys_sdEL(rL2ProjectionConfig%dabsError,10)))
        call output_line ("Iteraions : "//trim(sys_siL(rL2ProjectionConfig%iiterations,10)))
      end if
          
      ! Release the mass matrix
      call lsyssc_releaseMatrix (rmatrixMass)
    end do

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_writeucd (rcmdStatus,Sargs)
  
  !<description>
    ! Command: WRITEUCD.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=SYS_STRLEN) :: sfilename,smesh,stype,stoken,svecname,sname
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_vectorBlock), pointer :: p_rvectorBlock
    integer, dimension(:), pointer :: p_IdxVec
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
    type(t_meshHierarchy), pointer :: p_rmeshHierarchy
    type(t_ucdExport) :: rexport
    integer :: iparam
    logical :: bexists
    
    if (size(Sargs) .lt. 4) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Source vector and destination
    call cmdprs_getparam (Sargs,2,sfilename,.false.,.true.)
    call cmdprs_getparam (Sargs,3,stoken,.true.,.false.)
    call cmdprs_getparam (Sargs,4,stype,.true.,.false.)

    ! Get the mesh. Probably an index in a mesh hierarchy.
    call cmdprs_parseIndexString1D (stoken,smesh,p_IdxVec)
    if (.not. associated (p_IdxVec)) then
      p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,smesh,bexists=bexists)
      if (.not. bexists) then
        call output_line ("Unknown variable!")
        return
      end if
    else
      ! A mesh in a hierarchy.
      p_rmeshHierarchy => collct_getvalue_mshh(rcmdStatus%rcollection,smesh,bexists=bexists)
      if (.not. bexists) then
        deallocate(p_IdxVec)
        call output_line ("Unknown variable!")
        return
      end if
      
      if ((p_IdxVec(1) .lt. 1) .or. (p_IdxVec(1) .gt. p_rmeshHierarchy%nlevels)) then
        deallocate(p_IdxVec)
        call output_line ("Invalid mesh level!")
        return
      end if
      
      ! Get the mesh
      p_rtriangulation => p_rmeshHierarchy%p_Rtriangulations(p_IdxVec(1))
      deallocate(p_IdxVec)
    end if
    
    ! Open the UCD file.
    if (stype .eq. "VTK") then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)
    else if (stype .eq. "GMV") then
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)
    else
      call output_line ("Unknown output format!")
      return
    end if
    
    ! Loop through the parameters
    iparam = 5
    do while (iparam .le. size(Sargs))
      ! Next token
      call cmdprs_getparam (Sargs,iparam,stoken,.true.,.false.)
      iparam = iparam + 1
      
      if (stoken .eq. "--POINTDATASCALAR") then 
        
        ! Write scalar subvector for points
        
        if (iparam .ge. size(Sargs)) then
          ! Exit the loop, write the file.
          call output_line ("Invalid argument!")
          exit
        else
          ! Next tokens
          call cmdprs_getparam (Sargs,iparam,sname,.false.,.true.)
          call cmdprs_getparam (Sargs,iparam+1,stoken,.true.,.false.)
          iparam = iparam + 2
          
          ! Get the vector. Array specifier must be given.
          call cmdprs_parseIndexString1D (stoken,svecname,p_IdxVec)
          if ((svecname .eq. "") .or. (.not. associated(p_IdxVec))) then
            call output_line ("Invalid data vector!")
            exit
          end if
          
          p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svecname)
          if (.not. associated(p_rvectorBlock)) then
            deallocate(p_IdxVec)
            call output_line ("Invalid data vector!")
            exit
          end if
          
          if ((p_IdxVec(1) .lt. 1) .or. (p_IdxVec(1) .gt. p_rvectorBlock%nblocks)) then
            deallocate(p_IdxVec)
            call output_line ("Invalid component!")
            exit
          end if
          
          ! Write the block
          call spdp_projectToVertices (p_rvectorBlock%RvectorBlock(p_IdxVec(1)),p_Ddata1)
          call ucd_addVariableVertexBased (rexport,trim(sname),UCD_VAR_STANDARD,p_Ddata1)
          deallocate(p_Ddata1)
          
          deallocate(p_IdxVec)
          
        end if        

      else if (stoken .eq. "--POINTDATAVEC") then 
      
        ! Write vector field for points

        if (iparam .ge. size(Sargs)) then
          ! Exit the loop, write the file.
          call output_line ("Invalid argument!")
          exit
        else
          ! Next tokens
          call cmdprs_getparam (Sargs,iparam,sname,.false.,.true.)
          call cmdprs_getparam (Sargs,iparam+1,stoken,.true.,.false.)
          iparam = iparam + 2
          
          ! Get the vector. Array specifier must be given.
          call cmdprs_parseIndexString1D (stoken,svecname,p_IdxVec)
          if ((svecname .eq. "") .or. (.not. associated(p_IdxVec))) then
            call output_line ("Invalid data vector!")
            exit
          end if
          
          p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svecname)
          if (.not. associated(p_rvectorBlock)) then
            deallocate(p_IdxVec)
            call output_line ("Invalid data vector!")
            exit
          end if
          
          ! Up to three components
          nullify(p_Ddata1)
          nullify(p_Ddata2)
          nullify(p_Ddata3)
          
          if (size(p_IdxVec) .ge. 1) then
            if ((p_IdxVec(1) .lt. 1) .or. (p_IdxVec(1) .gt. p_rvectorBlock%nblocks)) then
              deallocate(p_IdxVec)
              call output_line ("Invalid component!")
              exit
            end if
            call spdp_projectToVertices (p_rvectorBlock%RvectorBlock(p_IdxVec(1)),p_Ddata1)
          end if  

          if (size(p_IdxVec) .ge. 2) then
            if ((p_IdxVec(2) .lt. 1) .or. (p_IdxVec(2) .gt. p_rvectorBlock%nblocks)) then
              deallocate(p_Ddata1)
              deallocate(p_IdxVec)
              call output_line ("Invalid component!")
              exit
            end if
            call spdp_projectToVertices (p_rvectorBlock%RvectorBlock(p_IdxVec(1)),p_Ddata2)
          end if  

          if (size(p_IdxVec) .ge. 3) then
            if ((p_IdxVec(2) .lt. 1) .or. (p_IdxVec(2) .gt. p_rvectorBlock%nblocks)) then
              deallocate(p_Ddata2)
              deallocate(p_Ddata1)
              deallocate(p_IdxVec)
              call output_line ("Invalid component!")
              exit
            end if
            call spdp_projectToVertices (p_rvectorBlock%RvectorBlock(p_IdxVec(1)),p_Ddata3)
          end if  
          
          ! Write the block
          select case (size(p_IdxVec))
          case (1)
            call ucd_addVarVertBasedVec (rexport,trim(sname),p_Ddata1)
            deallocate(p_Ddata1)
          case (2)
            call ucd_addVarVertBasedVec (rexport,trim(sname),p_Ddata1,p_Ddata2)
            deallocate(p_Ddata2)
            deallocate(p_Ddata1)
          case (3:)
            call ucd_addVarVertBasedVec (rexport,trim(sname),p_Ddata1,p_Ddata2,p_Ddata3)
            deallocate(p_Ddata3)
            deallocate(p_Ddata2)
            deallocate(p_Ddata1)
          end select
          
          deallocate(p_IdxVec)
          
        end if

      else if (stoken .eq. "--CELLDATASCALAR") then 

        ! Write scalar subvector for cells

        if (iparam .ge. size(Sargs)) then
          ! Exit the loop, write the file.
          call output_line ("Invalid argument!")
          exit
        else
          ! Next tokens
          call cmdprs_getparam (Sargs,iparam,sname,.false.,.true.)
          call cmdprs_getparam (Sargs,iparam+1,stoken,.true.,.false.)
          iparam = iparam + 2
          
          ! Get the vector. Array specifier must be given.
          call cmdprs_parseIndexString1D (stoken,svecname,p_IdxVec)
          if ((svecname .eq. "") .or. (.not. associated(p_IdxVec))) then
            call output_line ("Invalid data vector!")
            exit
          end if
          
          p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svecname)
          if (.not. associated(p_rvectorBlock)) then
            deallocate(p_IdxVec)
            call output_line ("Invalid data vector!")
            exit
          end if
          
          if ((p_IdxVec(1) .lt. 1) .or. (p_IdxVec(1) .gt. p_rvectorBlock%nblocks)) then
            deallocate(p_IdxVec)
            call output_line ("Invalid component!")
            exit
          end if
          
          ! Write the block
          call spdp_projectToCells (p_rvectorBlock%RvectorBlock(p_IdxVec(1)),p_Ddata1)
          call ucd_addVariableElementBased (rexport,trim(sname),UCD_VAR_STANDARD,p_Ddata1)
          deallocate(p_Ddata1)
          
          deallocate(p_IdxVec)
          
        end if        

      else
        call output_line ("Invalid parameter!")
        exit
      end if
    end do
    
    ! Write the file, done.
    call ucd_write(rexport)
    call ucd_release(rexport)
    
  end subroutine

end module
