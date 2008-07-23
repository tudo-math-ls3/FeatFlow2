!##############################################################################
!# ****************************************************************************
!# <name> prolrest </name>
!# ****************************************************************************
!#
!# <purpose>
!# Todo
!# </purpose>
!##############################################################################

PROGRAM prolrest
   
  USE prolrest2d_test1
  USE prolrest2d_test2
  USE prolrest2d_test3
  
  IMPLICIT NONE
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  CALL system_init()
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  CALL output_init ('./log/output.txt')

  ! The very second thing in every program: 
  ! Initialise the FEAT 2.0 storage management: 
  CALL storage_init(999, 100)

  ! Call 2D prol-rest test
  CALL prolrest2d_1()
  CALL prolrest2d_2()
  CALL prolrest2d_3() 

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk ()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
