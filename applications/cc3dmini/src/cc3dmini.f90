!##############################################################################
!# ****************************************************************************
!# <name> cc3dmini </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising a simple
!# stationary Navier-Stokes equation
!#
!#              $$- \nu Laplace(u) + u*grad(u) + \Nabla p = f $$
!#              $$ \Nable \cdot p = 0$$
!#
!# on a 3D domain for a 3D function $u=(u_1,u_2,u_3)$ and a pressure $p$.
!#
!# A linear block system with 4 solution components is set up,
!# discretised and solved.
!#
!# </purpose>
!##############################################################################

program cc3dmini

  use cc3dmini_method1
  
  implicit none
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
  call system_init()

  ! The very second thing in every program:
  ! Initialise the storage management:
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve. cc2d method 1:
  print *
  print *,'Calculating cc3dmini-Problem with method 1'
  print *,'------------------------------------------'
  call cc3dmini1

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  print *
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
