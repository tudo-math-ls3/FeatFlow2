!##############################################################################
!# ****************************************************************************
!# <name> prolrest </name>
!# ****************************************************************************
!#
!# <purpose>
!# Todo
!# </purpose>
!##############################################################################

program prolrest
   
  use prolrest2d_test1
  use prolrest2d_test2
  use prolrest2d_test3
  use prolrest2d_test4
  use prolrest2d_test5
  
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

  ! Call 2D prol-rest test
  !CALL prolrest2d_1()
  !CALL prolrest2d_2()
  !CALL prolrest2d_3()
  !CALL prolrest2d_4()
  call prolrest2d_5()

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
