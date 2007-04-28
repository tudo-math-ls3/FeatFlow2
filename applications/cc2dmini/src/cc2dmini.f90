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

PROGRAM cc2dmini

  USE cc2dmini_method1
  USE cc2dmini_method2
  
  IMPLICIT NONE
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()

  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  CALL storage_init(999, 100)

  ! Call the problem to solve. cc2d method 1:
  PRINT *
  PRINT *,'Calculating cc2dmini-Problem with method 1'
  PRINT *,'------------------------------------------'
  CALL cc2dmini1

  ! Call the problem to solve. cc2d method 2:
  PRINT *
  PRINT *,'Calculating cc2dmini-Problem with method 2'
  PRINT *,'------------------------------------------'
  CALL cc2dmini2

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  PRINT *
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
