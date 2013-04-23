!##############################################################################
!# Tutorial 001c: String manipulation routines
!##############################################################################

module tutorial001c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  implicit none
  private
  
  public :: start_tutorial001c

contains

  ! ***************************************************************************

  subroutine start_tutorial001c
  
    character(LEN=SYS_STRLEN) :: sstring

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 001c.")
    call output_separator (OU_SEP_MINUS)

    call output_line ("String mainipulation routines.")
    call output_line ("------------------------------")
    
    ! =================================
    ! Conversion to upper case / lower case
    ! =================================
    call output_line ("'This is a test' => " // sys_upcase ("'This is a test'") )
    call output_line ("'This is a test' => " // sys_lowcase("'This is a test'") )
    
    ! =================================
    ! Character replacement
    ! =================================
    call output_line ("'This is a test' => " // sys_charreplace("'This is a test'"," ","-") )

    ! =================================
    ! Dequoting. Partial output without line break
    ! =================================
    sstring = "'This is a test'"
    call output_line ( trim(sstring) // " =>", bnolinebreak=.true. )
    call sys_dequote ( sstring )
    call output_line ( " " // trim(sstring) )
    
    call output_lbrk()

  end subroutine

end module
