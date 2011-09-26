! A parser for standard IO commands: printf, sprintf.

module stdinoutparser

  use fsystem
  use genoutput
  use storage
  use collection
  use typedsymbol
  
  use commandparserbase

  implicit none
  
  private
  
  public :: sioprs_qualifyString
  
contains

  ! ***************************************************************************

  !<subroutine>

  subroutine sioprs_qualifyString (sstring,Rvariables)

  !<description>
    ! Resolves a string. Replaces all qualifiers in the string by the
    ! values in the variables.
  !</description>

  !<inputoutput>
    ! The string to resolve. Will be replaced by the fully qualified string.
    ! The string should not be dequoted.
    character(len=*), intent(inout) :: sstring
  !</inputoutput>
  
  !<input>
    ! List of variable values to be inserted into the string.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvariables
  !</input>
  
!</subroutine>
  
    ! local variables
    character(len=len(sstring)) :: soutput
    character(len=SYS_STRLEN) :: sformat
    integer :: iout, iin, ilen, istrlen, inextvar
    logical :: bescape,berror
    
    ! Loop through the input characters
    bescape = .false.
    iout = 1
    iin = 1
    inextvar = 1
    soutput = ""
    istrlen = len_trim(sstring)
    do 
      if (iin .gt. istrlen) exit
      
      if (bescape) then
        ! Take the next character as it is.
        soutput(iout:iout) = sstring(iin:iin)
        iout = iout + 1
      else
        ! What is the character?
        if (sstring(iin:iin) .eq. "\") then
          ! Escape the next character
          bescape = .true.
          
          ! Transfer the escape character, so later dequoting works correctly.
          soutput(iout:iout) = sstring(iin:iin)
          iout = iout + 1
          
        else if (sstring(iin:iin) .eq. "%") then
        
          ! Get the full format qualifier.
          call getFormat (sstring(iin:),sformat,ilen)
          
          ! Skip it in the input
          iin = iin + ilen - 1
          
          ! Plug in the next variable (correctly formatted) and continue
          ! in the output.
          if (present(Rvariables)) then
            if (inextvar .le. ubound(Rvariables,1)) then
              if ((Rvariables(inextvar)%ctype .ne. STYPE_INTEGER) .and. &
                  (Rvariables(inextvar)%ctype .ne. STYPE_DOUBLE) .and. &
                  (Rvariables(inextvar)%ctype .ne. STYPE_STRING)) then
                call output_line ("Ignoring invalid variable.");
                berror = .true.
              else
                call formatVariable (Rvariables(inextvar),soutput(iout:),sformat,ilen,berror)
                inextvar = inextvar + 1
              end if
              
              ! Ignore wrong qualifiers.
              if (.not. berror) then
                iout = iout + ilen
              end if
            end if
          end if
        else
          ! Transfer the character
          soutput(iout:iout) = sstring(iin:iin)
          iout = iout + 1
        end if
          
      end if
      
      bescape = .false.
      iin = iin + 1
    end do
    
    ! Return the string
    sstring = soutput
  
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine getFormat (sstring,sformat,ilen)

  !<description>
    ! Local routine. Extracts a format qualifier from a string.
    ! The qualifier must have the format "%... " ending with a space or
    ! the string end.
  !</description>

  !<input>
    ! String to analyse
    character(len=*), intent(in) :: sstring
  !</input>
  
  !<output>
    ! Format string.
    character(len=*), intent(out) :: sformat
    
    ! Length of the format string
    integer, intent(out) :: ilen
  !</output>
  
!</subroutine>

    ! Position of the first format specifier
    ilen = scan(sstring,"dfs")
    if (ilen .eq. 0) ilen = len_trim(sstring)
    
    sformat = sstring(1:ilen)
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine formatVariable (rvariable,sstring,sformat,ilen,berror)

  !<description>
    ! Local routine. Formats a variable according to the format string. 
  !</description>

  !<input>
    ! Variable to format.
    type(t_symbolValue), intent(in) :: rvariable
    
    ! Format string.
    character(len=*), intent(in) :: sformat
  !</input>
  
  !<output>
    ! Output string
    character(len=*), intent(inout) :: sstring
    
    ! Length of the output string
    integer, intent(out) :: ilen
    
    ! Set to TRUE if the format is wrong.
    logical, intent(out) :: berror
  !</output>
  
!</subroutine>

    integer :: ilenformat,ipostype,ios
    character (len=SYS_NAMELEN) :: snewformat
    character (len=SYS_STRLEN) :: stemp
    
    berror = .false.
    ilenformat = len_trim(sformat)
    
    ! Find the type
    ipostype = scan(sformat,"ldfs")
    
    ! Verify that the format is ok.
    if (verify(sformat(2:ipostype-1),"0123456789.") .ne. 0) then
      ! Incorrect string format
      berror = .true.
      return
    end if
    
    if (ipostype .gt. 2) then
      ! Get the format string
      snewformat = sformat(2:ipostype-1)
    else
      snewformat = ""
    end if

    ! What type should we write?
    if (sformat(ilenformat:ilenformat) .eq. "s") then
      ! Just a string.
      if (snewformat .eq. "") then
        ! Without length qualifier
        call cmdprs_dequoteStd (rvariable%svalue,sstring,ilen)
      else
        ! With length qualifier
        if (index(snewformat,".") .eq. 0) then
          read (snewformat,*) ilen
          snewformat = "(A"//trim(snewformat)//")"
          call cmdprs_dequoteStd (rvariable%svalue,stemp)
          write (sstring,snewformat,IOSTAT=ios) stemp
        else
          ! Otherwise: Format invalid.
          ilen = 0
          berror = .true.
        end if
      end if
      return
    end if
    
    if (sformat(ilenformat:ilenformat) .eq. "d") then
      ! Integer value.
      if (snewformat .ne. "") then
        snewformat = "(I"//trim(snewformat)//")"
        write (sstring,snewformat,IOSTAT=ios) rvariable%ivalue
      else
        write (sstring,*) rvariable%ivalue
        sstring = adjustl(sstring)
      end if
      ilen = len_trim(sstring)
      return
    end if

    if (sformat(ilenformat:ilenformat) .eq. "f") then
      ! Float value.
      if (snewformat .ne. "") then
        snewformat = "(F"//trim(snewformat)//")"
        write (sstring,snewformat,IOSTAT=ios) rvariable%dvalue
      else
        write (sstring,*) rvariable%dvalue
        sstring = adjustl(sstring)
      end if
      ilen = len_trim(sstring)
      return
    end if

    if (sformat(ilenformat:ilenformat) .eq. "e") then
      ! Float value.
      if (snewformat .ne. "") then
        snewformat = "(E"//trim(snewformat)//")"
        write (sstring,snewformat,IOSTAT=ios) rvariable%dvalue
      else
        write (sstring,*) rvariable%dvalue
        sstring = adjustl(sstring)
      end if
      ilen = len_trim(sstring)
      return
    end if

  end subroutine

end module
