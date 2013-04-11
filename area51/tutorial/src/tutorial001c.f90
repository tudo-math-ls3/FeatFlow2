!##############################################################################
!# Tutorial 001c: String manipulation routines
!##############################################################################

module tutorial001c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  implicit none
  private
  
  public :: start_tutorial001c

contains

  ! ***************************************************************************

  subroutine start_tutorial001c
  
    character(LEN=SYS_STRLEN) :: sstring
    integer :: ntokens

    ! Initialisation of the feat library, output system and memory management
    call system_init()
    call output_init ("")
    call storage_init(999, 100)

    ! Print a message
    call output_lbrk()
    call output_line ("Hello world. This is FEAT-2. Tutorial 001c.")
    call output_separator (OU_SEP_MINUS)

    call output_line ("String mainipulation routines.")
    call output_line ("------------------------------")
    
    ! Conversion to upper case / lower case
    call output_line ("'This is a test' => " // sys_upcase ("'This is a test'") )
    call output_line ("'This is a test' => " // sys_lowcase("'This is a test'") )
    
    ! Character replacement
    call output_line ("'This is a test' => " // sys_charreplace("'This is a test'"," ","-") )

    ! Dequoting. Partial output without line break
    sstring = "'This is a test'"
    call output_line ( trim(sstring) // " =>", bnolinebreak=.true. )
    call sys_dequote (sstring)
    call output_line ( " " // trim(sstring) )
    
    call output_lbrk()

    ! Clean up
    call storage_done()
    call output_done()
  
  end subroutine

end module
