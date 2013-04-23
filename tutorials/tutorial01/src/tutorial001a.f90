!##############################################################################
!# Tutorial 001a: Hello world
!##############################################################################

module tutorial001a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  implicit none
  private
  
  public :: start_tutorial001a

contains

  ! ***************************************************************************

  subroutine start_tutorial001a

    ! Print a message
    call output_line ("Hello world. This is FEAT-2. Tutorial 001a.")

  end subroutine

end module
