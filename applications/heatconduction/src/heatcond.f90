!##############################################################################
!# ****************************************************************************
!# <name> heatcond </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#   d/dt u  - div ( alpha*grad(u) )  +  beta*grad(u)  +  gamma*u  =  f
!#
!# on a 2D domain for a scalar function u=u(x,t).
!#
!# The first example (module heatcond_method5) is a very simple one-routine
!# solver for the head condition equation without any specials. The equation
!# is discretised with explicit Euler for a fixed number of time steps.
!# Boundary conditions and RHS are constant in time.
!#
!# The second example (module heatcond_method5) separates the different steps
!# of the solution process into different subroutines.
!# The communication is done using a problem-related structure. For the
!# communication with callback routines during the assembly, a
!# collection structure is set up. As solver, a simple Multigrid solver is
!# used.
!#
!# The equation itself is set up using parameters of a .DAT file
!# named 'data/heatcond.dat'. Here all parameters (alpha, beta, gamma, level,...)
!# can be found and changed by the user without recompiling the program.
!# </purpose>
!##############################################################################

PROGRAM heatcond

  USE heatcond_method1
  USE heatcond_method5
  
  IMPLICIT NONE

  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()

  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  CALL storage_init(999, 100)

  ! Call the problem to solve - method 5 = multigrid
  CALL output_lbrk()
  CALL output_line('Calculating heatcond-Problem with method 1')
  CALL output_line('------------------------------------------')
  CALL heatcond1
  
  ! Call the problem to solve - method 5 = multigrid
  CALL output_lbrk()
  CALL output_line('Calculating heatcond-Problem with method 5')
  CALL output_line('------------------------------------------')
  CALL heatcond5

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
