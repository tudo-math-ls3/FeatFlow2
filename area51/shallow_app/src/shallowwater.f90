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

program shallowwater
   
  use shallowwater2d
  
  implicit none
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
  call sys_init()
  sys_haltmode = SYS_HALT_THROWFPE
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  call output_init ('./log/output.txt')

  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)
 


  ! Call the problem to solve.
  call output_lbrk ()
  call output_line ('Calculating 2D shallow water problem')
  call output_line ('------------------------------------')
  call output_lbrk ()
  call shallowwater2d_0
  
 

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
  pause
  
end program
