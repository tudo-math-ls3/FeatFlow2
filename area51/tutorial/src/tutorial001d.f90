!##############################################################################
!# Tutorial 001d: Tokenizer
!##############################################################################

module tutorial001d

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  implicit none
  private
  
  public :: start_tutorial001d

contains

  ! ***************************************************************************

  subroutine start_tutorial001d
  
    character(LEN=SYS_STRLEN) :: sstring, stoken
    integer :: ntokens, itoken, istart

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("Hello world. This is FEAT-2. Tutorial 001d.")
    call output_separator (OU_SEP_MINUS)

    ! Count the number of tokens and print them.
    sstring = "0.1 0.2 'This is a test' 0.3"
    call output_line ("We now tokenize the string: """ // trim(sstring) // """.")
    call output_lbrk()
    
    call sys_countTokens (sstring,ntokens," ")
    call output_line ("Number of tokens: " // trim(sys_siL(ntokens,10)) )
    
    istart = 1
    itoken = 0
    
    do while (istart .ne. 0)
      itoken = itoken + 1
      call sys_getNextToken (sstring,stoken,istart)
      call output_line ("Token " // trim(sys_siL(itoken,10)) // ": " // trim(stoken) )
    end do

    call output_lbrk()

  end subroutine

end module
