!##############################################################################
!# ****************************************************************************
!# <name> bouss2dmini </name>
!# ****************************************************************************
!#
!# <purpose>
!# Todo
!# </purpose>
!##############################################################################

program bouss2dmini

  use bouss2dmini_mit_1
  use bouss2dmini_mit_2
  
  implicit none
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
  call sys_init()
  call output_init()

  ! The very second thing in every program:
  ! Initialise the storage management:
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve. 2d Boussinesq 1:
  call output_lbrk()
  call output_line('Calculating 2D Boussinesq-Problem - MIT - 1')
  call output_line('-------------------------------------------')
  call b2dm_mit_1

  ! Call the problem to solve. 2d Boussinesq 2:
  call output_lbrk()
  call output_line('Calculating 2D Boussinesq-Problem - MIT - 2')
  call output_line('-------------------------------------------')
  call b2dm_mit_2

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
