program featcmdline

  use fsystem
  use genoutput
  use storage
  use commandparserbase
  use commandparser

  ! Variables
  type(t_commandstatus) :: rcmdStatus

  ! The very first thing in every application:
  ! Initialise system-wide settings:

  call sys_init()

  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the converter
  call cmdprs_init(rcmdStatus)
  call cmdprs_parseterminal(rcmdStatus)
  call cmdprs_done(rcmdStatus)

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()

  ! Clean up the storage management, finish
  call storage_done()

end program
