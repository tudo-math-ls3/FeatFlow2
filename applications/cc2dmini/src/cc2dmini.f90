!##############################################################################
!# ****************************************************************************
!# <name> cc2dmini </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising a simple
!# stationary Navier-Stokes equation
!#
!#              $$- \nu Laplace(u) + u*grad(u) + \Nabla p = f $$
!#              $$ \Nable \cdot p = 0$$
!#
!# on a 2D domain for a 2D function $u=(u_1,u_2)$ and a pressure $p$.
!#
!# A linear block system with 3 solution components is set up,
!# discretised and solved.
!#
!# </purpose>
!##############################################################################

program cc2dmini

  use cc2dmini_method1
  use cc2dmini_method2
  
  implicit none
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
  call sys_init()

  ! The very second thing in every program:
  ! Initialise the storage management:
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve. cc2d method 1:
  print *
  print *,'Calculating cc2dmini-Problem with method 1'
  print *,'------------------------------------------'
  call cc2dmini1

  ! Call the problem to solve. cc2d method 2:
  print *
  print *,'Calculating cc2dmini-Problem with method 2'
  print *,'------------------------------------------'
  call cc2dmini2

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  print *
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
