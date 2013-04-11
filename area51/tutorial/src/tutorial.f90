!##############################################################################
!# Collection of tutorials.
!#
!# The main routine calls all tutorial applications.
!##############################################################################

program poisson

  ! Include the tutorial modules
  use tutorial001a
  use tutorial001b
  use tutorial001c
  use tutorial001d

  implicit none
  
  ! Call the tutorial modules.
  call start_tutorial001a
  call start_tutorial001b
  call start_tutorial001c
  call start_tutorial001d
  
end program
