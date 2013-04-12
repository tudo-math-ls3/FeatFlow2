!##############################################################################
!# Tutorial 002a: Read a mesh and write a VTK.
!##############################################################################

module tutorial002a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  implicit none
  private
  
  public :: start_tutorial002a

contains

  ! ***************************************************************************

  subroutine start_tutorial002a

    ! Print a message
    call output_line ("Hello world. This is FEAT-2. Tutorial 002a")

  end subroutine

end module
