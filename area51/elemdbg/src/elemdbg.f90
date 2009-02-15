!##############################################################################
!# ****************************************************************************
!# <name> elemdbg </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# </purpose>
!##############################################################################

program elemdbg
   
  use elemdbg1d_test1
  use elemdbg2d_test1
  use elemdbg3d_test1
  
  implicit none
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  call system_init()
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  call output_init ('./log/output.txt')

  ! The very second thing in every program: 
  ! Initialise the FEAT 2.0 storage management: 
  call storage_init(999, 100)
 
  ! 1D Element-Debugger
  call elemdbg1d_1()

  ! 2D Element-Debugger
  call elemdbg2d_1()

  ! 3D Element-Debugger
  call elemdbg3d_1()

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
