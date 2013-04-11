!##############################################################################
!# Tutorial 001a: Hello world
!##############################################################################

module tutorial001a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  implicit none

contains

  ! ***************************************************************************

  subroutine start_tutorial001a

    ! Initialisation of the feat library, output system and memory management
    call system_init()
    call output_init ("")
    call storage_init(999, 100)

    ! Print a message
    call output_line ("Hello world. This is FEAT-2.")

    ! Clean up
    call storage_done()
  
  end subroutine

end module
