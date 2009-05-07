!##############################################################################
!# ****************************************************************************
!# <name> shallowwater </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the 
!# shallowwater equations.
!#
!# </purpose>
!##############################################################################

PROGRAM shallowwater
   
  USE shallowwater2d
  
  IMPLICIT NONE
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()
  sys_haltmode = SYS_HALT_THROWFPE
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  CALL output_init ('./log/output.txt')

  ! The very second thing in every program: 
  ! Initialise the FEAT 2.0 storage management: 
  CALL storage_init(999, 100)
 


  ! Call the problem to solve.
  CALL output_lbrk ()
  CALL output_line ('Calculating 2D shallow water problem')
  CALL output_line ('------------------------------------')
  CALL output_lbrk ()
  CALL shallowwater2d_0
  
 

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk ()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
  PAUSE
  
END PROGRAM
