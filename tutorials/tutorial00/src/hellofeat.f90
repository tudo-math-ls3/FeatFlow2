!####################################################################
!# The simplest FEAT2 application. Hello world in FEAT2.
!####################################################################

program HelloFeat

  use fsystem
  use genoutput

  implicit none

  ! Initialisation
  call sys_init()
  call output_init("")

  ! Print a message
  call output_line ("Hello world. This is FEAT-2.")

  ! Cleanup
  call output_done()

end program
