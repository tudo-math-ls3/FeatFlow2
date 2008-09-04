!##############################################################################
!# ****************************************************************************
!# <name> bouss2dmini </name>
!# ****************************************************************************
!#
!# <purpose>
!# Todo
!# </purpose>
!##############################################################################

PROGRAM bouss2dmini

  USE bouss2dmini_mit_1
  USE bouss2dmini_mit_2
  
  IMPLICIT NONE
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()
  CALL output_init()

  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  CALL storage_init(999, 100)

  ! Call the problem to solve. 2d Boussinesq 1:
  CALL output_lbrk()
  CALL output_line('Calculating 2D Boussinesq-Problem - MIT - 1')
  CALL output_line('-------------------------------------------')
  CALL b2dm_mit_1

  ! Call the problem to solve. 2d Boussinesq 2:
  CALL output_lbrk()
  CALL output_line('Calculating 2D Boussinesq-Problem - MIT - 2')
  CALL output_line('-------------------------------------------')
  CALL b2dm_mit_2

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
