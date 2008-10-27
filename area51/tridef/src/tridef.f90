!##############################################################################
!# ****************************************************************************
!# <name> tridef </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#              - Laplace(u) = f
!#
!# for a scalar function u.
!#
!# There are a couple of examples provided how to solve this problem.
!# Each example has its own characteristics how to solve the problem.
!# 
!# </purpose>
!##############################################################################

program tridef
   
  USE tridef2d_method0
  
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
  ! Call the problem to solve. Poisson 1D method 1 - multigrid:

  ! Call the problem to solve. Poisson 2D method 1 - simple:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 0 - simple')
  CALL output_line ('-----------------------------------------------------')
  CALL tridef2d_simple
  
  
    ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
