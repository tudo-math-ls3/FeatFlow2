!##############################################################################
!# Tutorial 001b: Basic output commands
!##############################################################################

module tutorial001b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  implicit none
  private
  
  public :: start_tutorial001b

contains

  ! ***************************************************************************

  subroutine start_tutorial001b

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("Hello world. This is FEAT-2. Tutorial 001b.")    
    call output_separator (OU_SEP_MINUS)
    
    call output_line ("Let us print a table -- right-adjusted.")
    call output_lbrk()
    call output_line ("  n          x^n")
    call output_line ("-----------------")
    call output_line (sys_si(0,3) // "|" // sys_sd(2.0_DP ** 0,10) // "|")
    call output_line (sys_si(1,3) // "|" // sys_sd(2.0_DP ** 1,10) // "|")
    call output_line (sys_si(4,3) // "|" // sys_sd(2.0_DP ** 4,10) // "|")
    call output_lbrk()

    call output_separator (OU_SEP_MINUS)
    call output_line ("Let us print a table -- left-adjusted.")
    call output_lbrk()
    call output_line ("  n          x^n")
    call output_line ("-----------------")
    call output_line (sys_siL(0,3) // "|" // sys_sdL(2.0_DP ** 0,10) // "|")
    call output_line (sys_siL(1,3) // "|" // sys_sdL(2.0_DP ** 1,10) // "|")
    call output_line (sys_siL(4,3) // "|" // sys_sdL(2.0_DP ** 4,10) // "|")
    call output_lbrk()

    call output_separator (OU_SEP_MINUS)
    call output_line ("Let us print a table -- nicely adjusted.")
    call output_lbrk()
    call output_line ("  n              x^n")
    call output_line ("---------------------")
    call output_line (trim(sys_si(0,3)) // "|" // sys_adjustr(sys_sdL(2.0_DP ** 0,10),16) // "|" )
    call output_line (trim(sys_si(1,3)) // "|" // sys_adjustr(sys_sdL(2.0_DP ** 1,10),16) // "|" )
    call output_line (trim(sys_si(4,3)) // "|" // sys_adjustr(sys_sdL(2.0_DP ** 4,10),16) // "|" )
    call output_lbrk()

    call output_separator (OU_SEP_MINUS)
    call output_line ("Let us print a table -- scientific.")
    call output_lbrk()
    call output_line ("  n              x^n")
    call output_line ("---------------------")
    call output_line (trim(sys_si(0,3)) // "|" // trim(sys_sdEL(2.0_DP ** 0,10)) // "|" )
    call output_line (trim(sys_si(1,3)) // "|" // trim(sys_sdEL(2.0_DP ** 1,10)) // "|" )
    call output_line (trim(sys_si(4,3)) // "|" // trim(sys_sdEL(2.0_DP ** 4,10)) // "|" )
    call output_lbrk()

  end subroutine

end module
