!##############################################################################
!# ****************************************************************************
!# <name> cc2dmedium </name>
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

PROGRAM cc2dmedium

  USE cc2dmedium_method2
  
  IMPLICIT NONE
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()
  
  ! Initialise log file for output.
  CALL output_init ('log/output.txt')
  
  ! The very second thing in every program: 
  ! Initialise the storage management: 
  CALL storage_init(999, 100)
  
  ! Initialise the parser
  CALL fparser_init ()

  ! Call the problem to solve. 
  CALL output_lbrk ()
  CALL output_line ('Calculating cc2dmedium-Problem')
  CALL output_separator (OU_SEP_MINUS)
  
  CALL cc2dmedium2 ()

  ! Release the parser
  CALL fparser_done ()

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk ()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
