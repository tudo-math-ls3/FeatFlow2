module commandparserbase

  use fsystem
  use collection
  use ISO_VARYING_STRING

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
  
  public :: t_commandstatus
  
  public :: cmdprs_dequoteStd
  public :: cmdprs_dequote
  public :: cmdprs_getSection
  public :: cmdprs_getSymbolSection
  public :: cmdprs_complexSplit
  public :: cmdprs_nexttoken
  public :: cmdprs_counttokens
  public :: cmdprs_commandmatch
  public :: cmdprs_parseIndexString1D
  public :: cmdprs_removeTrail

contains

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

  !************************************************************************

  elemental subroutine cmdprs_dequoteStd (ssource,sdest,ilength)
  !<description>
    ! Remove all quotation marks from a string.
  !</description>
  
  !<input>
    ! A source string.
    character(len=*), intent(in) :: ssource
  !</input>
  
  !<output>
    ! Destination string. Non-escaped quotation marks are removed.
    character(len=*), intent(out) :: sdest
    
    ! OPTIONAL: Length of the string
    integer, intent(out), optional :: ilength
  !</output>
  
    character(len=len(ssource)) :: stemp
    integer :: i, idest, ilen
    character :: squotechar, scurrentchar
    logical :: bescape
    
    ! Copy the string, remove leading/trailing spaces.
    stemp = trim(ssource)
    
    ! Initialise destination
    sdest = ""
    
    ! Loop through the string. Copy characters.
    ! Replace quotation marks.
    i = 1
    idest = 0
    bescape = .false.
    squotechar = ' '
    ilen = len_trim(stemp)
    do 
      if (i .gt. ilen) exit
      
      ! Get the current character.
      scurrentchar = stemp(i:i)
      
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
      sdest(idest:idest) = scurrentchar

      i = i+1
    end do
    
    ! From the destination, take the actual substring.
    if (present(ilength)) then
      ilength = idest
    end if
    
  end subroutine 

  ! ***************************************************************************

  !<subroutine>

  subroutine cmdprs_getSymbolSection (rcollection,ssymbol,inestlevel,ssection,bfound)

  !<description>
    ! Tries to find the name of the section containing a symbol.
  !</description>

  !<inputoutput>
    ! Collection containing symbols.
    type(t_collection), intent(in) :: rcollection
  !</inputoutput>
  
  !<input>
    ! Name of the symbol (variable,...)
    character(len=*), intent(in) :: ssymbol
    
    ! Level of nesting
    integer, intent(in) :: inestlevel
  !</input>
  
  !<output>
    ! Name of the section.
    character(len=*), intent(out) :: ssection
    
    ! OPTIONAL: Returns if the symbol was found or not.
    logical, intent(out), optional :: bfound
  !</output>
  
!</subroutine>
    
    ! local variables
    integer :: ilv

    if (present(bfound)) then
      bfound = .false.
    end if
    ssection = ""
    
    ! Start searching from the deepest nesting
    do ilv = inestlevel,0,-1
      call cmdprs_getSection (ilv,ssection)
      
      if (collct_queryvalue(rcollection,ssymbol,ssectionName=ssection) .gt. 0) then
        ! Found.
        if (present(bfound)) then
          bfound = .true.
        end if
        return
      end if
    end do
    
    ! Not found. Return the innermost section
    call cmdprs_getSection (inestlevel,ssection)
  
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine cmdprs_getSection (inestlevel,ssection)

  !<description>
    ! Returns the section corresponding to a nesting level.
  !</description>

  !<input>
    ! Level of nesting
    integer, intent(in) :: inestlevel
  !</input>
  
  !<output>
    ! Name of the section.
    character(len=*), intent(out) :: ssection
  !</output>
  
!</subroutine>
    
    ! Return the section.
    if (inestlevel .gt. 0) then
      ! local variables
      ssection = trim(sys_siL(inestlevel,10))
    else
      ! Unnamed section, global variables
      ssection = ""
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
    !    Group 2: "!+-*/%()~,;{}^<>="
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
          call trigger_chargroup2 (ccurrchargroup,"!+-*/%()~,;{}^<>=",btriggered,ssource,ipos,&
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

  subroutine cmdprs_nexttoken (ssource,istart,iend,ilength,bback,cseparator)
  
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
    
    ! OPTIONAL: Character which is used to separate the tokens.
    ! Default s char(0).
    character, intent(in), optional :: cseparator
  !</inputoutput>
  
  !</subroutine>
  
    integer :: ilen
    logical :: back
    character :: csep
    
    if (present(ilength)) then
      ilen = ilength
    else
      ilen = len_trim(ssource)
    end if
    
    csep = char(0)
    if (present(cseparator)) csep = cseparator
    
    back = .false.
    if (present(bback)) back = bback
    
    if (.not. back) then
      
      ! Forward mode
      if (istart .eq. 0) then
        
        ! Find first token
        istart = 1
        iend = index (ssource,csep)
        if (iend .eq. 0) then
          iend = ilen
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
          iend = index (ssource(istart:),csep)
          if (iend .eq. 0) then
            iend = ilen
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
        istart = index (ssource,csep,.true.) + 1
        
      else
      
        ! Find next token
        iend = istart - 2
        if (iend .lt. 1) then
          ! Finish
          istart = 0
          iend = 0
        else
          ! Find next end.
          istart = index (ssource(1:iend),csep,.true.) + 1
        end if
      end if
    
    end if

  end subroutine

  ! ***************************************************************************

  !</subroutine>

  integer function cmdprs_counttokens (ssource,cseparator)
  
  !<description>
    ! Calculates the number of tokens in ssource.
  !</description>
  
  !<input>
    ! A source string. Must have been split with cmdprs_complexSplit.
    character(len=*), intent(in) :: ssource

    ! OPTIONAL: Character which is used to separate the tokens.
    ! Default s char(0).
    character, intent(in), optional :: cseparator
  !</input>
  
  !<result>
    ! Number of tokens. =0 if ssource does not contain text.
  !</result>
  
  !</subroutine>
  
    integer :: i,ntokens
    character :: csep
    
    csep = char(0)
    if (present(cseparator)) csep = cseparator

    cmdprs_counttokens = 0
    
    ! Return if there is no text.
    if (ssource .eq. "") return
  
    ! Count the number of \0 characters
    ntokens = 0
    do i = 1,len_trim(ssource)
      if (ssource(i:i) .eq. csep) then
        ntokens = ntokens + 1
      end if
    end do
    
    ! Add 1 for the last token.
    ntokens = ntokens + 1
    
    cmdprs_counttokens = ntokens
    
  end function

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
        !else if (ssource(icurrenttokenstart:icurrenttokenend) .eq. ";") then
          ! Done. ";" at the end is ignored.
        !  return
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
      
      !if (ssource(icurrenttokenstart:icurrenttokenend) .eq. ";") then
      !  ! Ignore ";" at the end.
      !  cycle
      !end if
      
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

  subroutine cmdprs_removeTrail (sline,btrailexists)
  
  !<description>
    ! Removes a trailing semicolon from a line.
  !</description>
  
  !<inputoutput>
    ! Splitted line of symbols. If there is a trailing semicolon, it will be
    ! removed.
    character(len=*), intent(inout) :: sline
    
    ! Returns whether or not there was a trailing semicolon.
    logical, intent(out) :: btrailexists
  !</inputoutput>
  
  !</subroutine>
  
    integer :: istart,iend
    
    istart=0
    iend=0
    
    ! Take a lok at the last token.
    call cmdprs_nexttoken (sline,istart,iend,len(sline),.true.)
    
    ! Semicolon?
    btrailexists = sline(istart:iend) .eq. ";"
    if (btrailexists) then
      ! Remove that trail.
      if (istart .gt. 1) then
        sline(istart:) = ""
      else
        sline = ""
      end if
    end if
  
  end subroutine

end module
