!####################################################################
!# Collection of tutorials.
!#
!# The main routine calls all tutorial applications.
!####################################################################

program tutorial

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  ! Include the tutorial modules
  use tutorial001a
  use tutorial001b
  use tutorial001c
  use tutorial001d

  implicit none
  
  ! =================================================================
  ! Main program
  ! =================================================================
  
  ! -----------------------------------------------------------------
  ! Initialisation of the feat library, output system and 
  ! memory management
  ! -----------------------------------------------------------------
  call system_init()
  call output_init ("")
  call storage_init(999, 100)

  ! -----------------------------------------------------------------
  ! Call the tutorial modules.
  ! -----------------------------------------------------------------
  call start_tutorial001a
  call start_tutorial001b
  call start_tutorial001c
  call start_tutorial001d
  
  ! -----------------------------------------------------------------
  ! Clean up
  ! -----------------------------------------------------------------
  call storage_done()
  call output_done()

end program
