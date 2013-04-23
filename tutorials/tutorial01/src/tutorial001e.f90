!##############################################################################
!# Tutorial 001e: Environment variables
!##############################################################################

module tutorial001e

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  implicit none
  private
  
  public :: start_tutorial001e

contains

  ! ***************************************************************************

  subroutine start_tutorial001e
  
    character(LEN=2048) :: sresult

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 001e.")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Print the content of the PATH
    ! variable.
    ! =================================
    if (.not. sys_getenv_string("PATH", sresult)) then
      call output_line ("Cannot get the PATH variable..")
    else
      call output_line ("$PATH = " // trim(sresult))
    end if
    
  end subroutine

end module
